"""Container-side preparation and execution of generated MHD cases."""

from __future__ import annotations

import json
import re
import shutil
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
        dry_run: bool = False,
    ) -> Path:
        """Run a previously prepared coupled OpenFOAM/Elmer case."""
        case_directory = self.case_directory(self._validate_case_name(name))
        marker = self._read_marker(case_directory)
        ranks = marker.get("ranks")
        self._validate_ranks(ranks)
        self._validate_prepared_case(case_directory, ranks)

        commands: List[Sequence[str]] = [
            (
                "mpirun",
                "-n",
                str(ranks),
                "mdhLinearHall",
                "-parallel",
                ":",
                "-n",
                str(ranks),
                "ElmerSolver_mpi",
                "case.sif",
            )
        ]
        if postprocess:
            commands.extend(
                (
                    ("reconstructPar", "-case", "."),
                    ("foamToVTK", "-case", "."),
                )
            )

        if dry_run:
            print(f"Would run {case_directory}")
            for command in commands:
                print(f"  run: {' '.join(command)}")
            return case_directory

        for command in commands:
            self._run(command, case_directory)
        return case_directory

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
            ("potentialFoam",),
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
