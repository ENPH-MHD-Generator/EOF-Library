"""Container-side preparation and execution of generated MHD cases."""

from __future__ import annotations

import json
import re
import shutil
import os
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from .config import CaseConfig, load_case_config
from .preparer import CasePreparer, PreparationError


CASE_NAME_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
SUBDOMAIN_PATTERN = re.compile(r"(numberOfSubdomains\s+)\d+(\s*;)")


class RuntimeCommandError(RuntimeError):
    """Raised when an external case-preparation or solver command fails."""


class MhdRuntime:
    """Prepare and execute cases rooted beneath a fixed runs directory."""

    MARKER_NAME = ".mhd-prepared.json"

    def __init__(
        self,
        base_case_directory: Path,
        runs_directory: Path = Path("/runs"),
        *,
        experiments_directory: Path = Path("/experiments"),
    ) -> None:
        self.base_case_directory = base_case_directory.resolve()
        self.runs_directory = runs_directory.resolve()
        self.experiments_directory = experiments_directory.resolve()
        self.case_preparer = CasePreparer(self.base_case_directory)

    def available_experiments(self) -> List[Path]:
        """Return every YAML experiment currently visible in the input mount."""
        if not self.experiments_directory.is_dir():
            return []
        return sorted(
            path
            for path in self.experiments_directory.rglob("*")
            if path.is_file()
            and not path.is_symlink()
            and path.suffix.lower() in {".yaml", ".yml"}
            and not any(
                part.startswith(".")
                for part in path.relative_to(self.experiments_directory).parts
            )
        )

    def available_cases(self) -> List[str]:
        """Return prepared case names that pass the normal run preflight checks."""
        if not self.runs_directory.is_dir():
            return []

        available: List[str] = []
        candidates = sorted(
            self.runs_directory.iterdir(), key=lambda path: path.name
        )
        for candidate in candidates:
            if not candidate.is_dir() or not CASE_NAME_PATTERN.fullmatch(candidate.name):
                continue
            try:
                case_directory = self.case_directory(candidate.name)
                marker = self._read_marker(case_directory)
                ranks = marker.get("ranks")
                self._validate_ranks(ranks)
                self._validate_prepared_case(case_directory, ranks)
            except (OSError, PreparationError):
                continue
            available.append(candidate.name)
        return available

    def prepare(
        self,
        config_path: Path,
        *,
        name: Optional[str] = None,
        generate_mesh: bool = False,
        ranks: int = 2,
        force: bool = False,
        dry_run: bool = False,
    ) -> Path:
        """Compile and convert a case into an execution-ready directory."""
        config_path = config_path.resolve()
        config = load_case_config(config_path)
        case_name = self._validate_case_name(name or config_path.stem)
        self._validate_ranks(ranks)
        case_directory = self.case_directory(case_name)

        if dry_run:
            mode = "generate a fresh mesh" if generate_mesh else "reuse the base mesh"
            print(f"Would prepare {case_directory}")
            print(f"  configuration: {config_path}")
            print(f"  mesh: {mode}")
            print(f"  MPI ranks per solver: {ranks}")
            for command in self.preparation_commands(ranks):
                print(f"  run: {' '.join(command)}")
            for item in self.case_preparer.describe(
                config,
                case_directory,
                generate_mesh=generate_mesh,
                reuse_mesh=not generate_mesh,
            ):
                print(f"  write: {item}")
            return case_directory

        if case_directory.exists() and not force:
            raise PreparationError(
                f"Case '{case_name}' already exists; pass --force to replace it"
            )
        self.runs_directory.mkdir(parents=True, exist_ok=True)
        workspace = Path(
            tempfile.mkdtemp(prefix=f".{case_name}.preparing-", dir=self.runs_directory)
        )
        staged_case = workspace / "case"
        try:
            self.case_preparer.prepare(
                config,
                staged_case,
                generate_mesh=generate_mesh,
                reuse_mesh=not generate_mesh,
            )
            self._set_openfoam_ranks(staged_case, ranks)
            for command in self.preparation_commands(ranks):
                self._run(command, staged_case, display_name=case_name)
            (staged_case / "ELMERSOLVER_STARTINFO").write_text(
                "case.sif\n", encoding="utf-8"
            )
            self._validate_prepared_case(staged_case, ranks)
            self._write_marker(staged_case, config_path, config, ranks)
            self._promote(staged_case, case_directory, workspace)
        finally:
            shutil.rmtree(workspace, ignore_errors=True)
        return case_directory

    def run(
        self,
        name: str,
        *,
        postprocess: bool = True,
        cores: Optional[int] = None,
        dry_run: bool = False,
    ) -> Path:
        """Run a previously prepared coupled OpenFOAM/Elmer case."""
        case_directory = self.case_directory(self._validate_case_name(name))
        marker = self._read_marker(case_directory)
        ranks = marker.get("ranks")
        self._validate_ranks(ranks)
        self._validate_prepared_case(case_directory, ranks)

        if cores is not None:
            available = os.cpu_count() or 1
            if isinstance(cores, bool) or not isinstance(cores, int) or not 1 <= cores <= available:
                raise PreparationError(f"--cores must be between 1 and {available}")
        # OpenFOAM and Elmer take turns: while one computes, the other waits in
        # the coupler's sleeping poll. Both solvers can therefore share the same
        # cores (N + N processes on N cores). Oversubscription must be allowed,
        # binding off, and MPI told to yield when idle so waiting ranks inside
        # MPI calls give up the CPU instead of busy-polling.
        pinning: List[str] = (
            ["taskset", "-c", f"0-{cores - 1}"] if cores is not None else []
        )
        solver_command = pinning + [
            "mpirun",
            "--oversubscribe",
            "--bind-to",
            "none",
            "--mca",
            "mpi_yield_when_idle",
            "1",
            "-n",
            str(ranks),
            "mdhLinearHall",
            "-parallel",
            ":",
            "-n",
            str(ranks),
            "ElmerSolver_mpi",
            "case.sif",
        ]
        workers = cores if cores is not None else ranks

        if dry_run:
            print(f"Would run {case_directory}")
            print(f"  run: {' '.join(solver_command)}")
            if postprocess:
                print(f"  run: reconstructPar and foamToVTK over the time directories, {workers} at a time")
            return case_directory

        self._run(solver_command, case_directory)
        if postprocess:
            self._postprocess(case_directory, workers, pinning)
        return case_directory

    def _postprocess(self, case_directory: Path, workers: int, pinning: Sequence[str]) -> None:
        """Reconstruct fields and write VTK, splitting time directories across workers.

        Each process handles its own subset of times, so the outputs are the same
        files a single serial reconstructPar / foamToVTK run writes.
        """
        processor_times = self._time_directories(case_directory / "processor0")
        self._parallel_over_times(
            "reconstructPar",
            [t for t in processor_times if float(t) != 0.0],
            workers,
            case_directory,
            pinning,
        )
        self._parallel_over_times(
            "foamToVTK", self._time_directories(case_directory), workers, case_directory, pinning
        )

    @staticmethod
    def _time_directories(directory: Path) -> List[str]:
        times = []
        for entry in directory.iterdir() if directory.is_dir() else []:
            if not entry.is_dir():
                continue
            try:
                value = float(entry.name)
            except ValueError:
                continue
            times.append((value, entry.name))
        return [name for _, name in sorted(times)]

    def _parallel_over_times(
        self,
        tool: str,
        times: Sequence[str],
        workers: int,
        case_directory: Path,
        pinning: Sequence[str],
    ) -> None:
        if not times:
            return
        if shutil.which(tool) is None:
            raise RuntimeCommandError(f"Required command '{tool}' is not available in the container")
        count = max(1, min(workers, len(times)))
        # Round-robin keeps the groups similar in size
        groups = [list(times[i::count]) for i in range(count)]
        name = case_directory.name
        print(f"[{name}] $ {tool} over {len(times)} time(s) in {count} parallel process(es)", flush=True)
        processes = []
        for index, group in enumerate(groups):
            log_path = case_directory / f"log.{tool}.{index}"
            command = list(pinning) + [tool, "-case", ".", "-time", ",".join(group)]
            log = open(log_path, "w")
            processes.append((subprocess.Popen(command, cwd=case_directory, stdout=log, stderr=subprocess.STDOUT), log, log_path, command))
        failures = []
        for process, log, log_path, command in processes:
            process.wait()
            log.close()
            if process.returncode != 0:
                failures.append((log_path, command, process.returncode))
        if failures:
            log_path, command, code = failures[0]
            tail = log_path.read_text(errors="replace").splitlines()[-15:]
            raise RuntimeCommandError(
                f"Command failed with exit code {code}: {' '.join(command)}\n"
                + "\n".join(tail)
                + f"\n(full log: {log_path})"
            )

    def case_directory(self, name: str) -> Path:
        """Resolve a safe case name beneath the immutable runs root."""
        destination = (self.runs_directory / name).resolve()
        try:
            destination.relative_to(self.runs_directory)
        except ValueError as exc:
            raise PreparationError("Case output must remain beneath /runs") from exc
        return destination

    @staticmethod
    def preparation_commands(ranks: int) -> List[Sequence[str]]:
        return [
            ("gmshToFoam", "channel.msh"),
            # Initial velocity from potential flow. -pName points at a field that
            # does not exist so the pressure BCs are inferred from U (the
            # absolute-pressure field p is not a kinematic pressure)
            ("potentialFoam", "-pName", "pPotential"),
            # potentialFoam also writes its volumetric flux; the compressible
            # solver would read it as the mass flux
            ("rm", "-f", "0/phi"),
            ("decomposePar",),
            (
                "ElmerGrid",
                "14",
                "2",
                "channel.msh",
                "-out",
                "meshElmer",
                "-autoclean",
                "-merge",
                "1e-8",
                "-removeunused",
            ),
            ("ElmerGrid", "2", "2", "meshElmer", "-metis", str(ranks)),
        ]

    @staticmethod
    def _run(
        command: Sequence[str],
        case_directory: Path,
        *,
        display_name: Optional[str] = None,
    ) -> None:
        executable = shutil.which(command[0])
        if executable is None:
            raise RuntimeCommandError(
                f"Required command '{command[0]}' is not available in the container"
            )
        print(f"[{display_name or case_directory.name}] $ {' '.join(command)}", flush=True)
        try:
            subprocess.run(command, cwd=case_directory, check=True)
        except subprocess.CalledProcessError as exc:
            raise RuntimeCommandError(
                f"Command failed with exit code {exc.returncode}: {' '.join(command)}"
            ) from exc

    @staticmethod
    def _validate_case_name(name: str) -> str:
        if not CASE_NAME_PATTERN.fullmatch(name) or name in {".", ".."}:
            raise PreparationError(
                "Case names may contain only letters, numbers, '.', '_', and '-', "
                "and must start with a letter or number"
            )
        return name

    @staticmethod
    def _validate_ranks(ranks: object) -> None:
        if isinstance(ranks, bool) or not isinstance(ranks, int) or ranks <= 0:
            raise PreparationError("MPI ranks must be a positive integer")

    @staticmethod
    def _set_openfoam_ranks(case_directory: Path, ranks: int) -> None:
        dictionary = case_directory / "system" / "decomposeParDict"
        content = dictionary.read_text(encoding="utf-8")
        updated, replacements = SUBDOMAIN_PATTERN.subn(
            rf"\g<1>{ranks}\g<2>", content, count=1
        )
        if replacements != 1:
            raise PreparationError(
                f"Could not update numberOfSubdomains in '{dictionary}'"
            )
        dictionary.write_text(updated, encoding="utf-8")

    def _write_marker(
        self,
        case_directory: Path,
        config_path: Path,
        config: CaseConfig,
        ranks: int,
    ) -> None:
        marker: Dict[str, object] = {
            "version": 1,
            "prepared_at": datetime.now(timezone.utc).isoformat(),
            "configuration": str(config_path),
            "ranks": ranks,
            "config": config.to_dict(),
        }
        (case_directory / self.MARKER_NAME).write_text(
            json.dumps(marker, indent=2) + "\n", encoding="utf-8"
        )

    @staticmethod
    def _promote(staged_case: Path, destination: Path, workspace: Path) -> None:
        """Atomically expose a fully prepared case, restoring the old one on error."""
        previous = workspace / "previous"
        if destination.exists():
            destination.replace(previous)
        try:
            staged_case.replace(destination)
        except Exception:
            if previous.exists():
                previous.replace(destination)
            raise

    def _read_marker(self, case_directory: Path) -> Dict[str, object]:
        marker_path = case_directory / self.MARKER_NAME
        if not marker_path.is_file():
            raise PreparationError(
                f"Case '{case_directory.name}' is not prepared; run 'mhd prepare' first"
            )
        try:
            data = json.loads(marker_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise PreparationError(f"Invalid preparation marker '{marker_path}'") from exc
        if not isinstance(data, dict) or data.get("version") != 1:
            raise PreparationError(f"Unsupported preparation marker '{marker_path}'")
        return data

    @staticmethod
    def _validate_prepared_case(case_directory: Path, ranks: int) -> None:
        required = [
            case_directory / "case.sif",
            case_directory / "ELMERSOLVER_STARTINFO",
            case_directory / "channel.msh",
            case_directory / "meshElmer",
            case_directory / "constant" / "polyMesh",
        ]
        required.extend(case_directory / f"processor{index}" for index in range(ranks))
        missing = [str(path) for path in required if not path.exists()]
        if missing:
            raise PreparationError(
                "Prepared case is missing required paths: " + ", ".join(missing)
            )
