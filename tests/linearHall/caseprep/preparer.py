"""High-level orchestration for producing complete, isolated case directories."""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

import yaml

from .boundaries import elmer_boundary_indices
from .config import CaseConfig, ConfigError, load_case_config
from .mesh import ProceduralMeshGenerator
from .renderers import ElmerCaseRenderer, OpenFoamCaseRenderer, render_coupling_properties


class PreparationError(RuntimeError):
    """Raised when a case cannot be materialized safely."""


class CasePreparer:
    """Compile a validated configuration into a self-contained solver case."""

    def __init__(
        self,
        base_case_directory: Path,
        *,
        elmer_renderer: Optional[ElmerCaseRenderer] = None,
        openfoam_renderer: Optional[OpenFoamCaseRenderer] = None,
        mesh_generator: Optional[ProceduralMeshGenerator] = None,
    ) -> None:
        self.base_case_directory = base_case_directory.resolve()
        self.elmer_renderer = elmer_renderer or ElmerCaseRenderer()
        self.openfoam_renderer = openfoam_renderer or OpenFoamCaseRenderer()
        self.mesh_generator = mesh_generator or ProceduralMeshGenerator()

    def describe(
        self,
        config: CaseConfig,
        destination: Path,
        *,
        generate_mesh: bool = False,
        reuse_mesh: bool = False,
    ) -> List[str]:
        self._validate_options(generate_mesh, reuse_mesh)
        self._validate_base_case(config, reuse_mesh=reuse_mesh)
        files = ["case.sif", "constant/couplingProperties", "resolved-config.yaml", "manifest.json"]
        files.extend(f"0/{name}" for name in self.openfoam_renderer.render(config))
        files.extend(("constant/", "system/"))
        if generate_mesh or reuse_mesh:
            files.append("channel.msh")
        return [f"{destination.resolve()} <- {item}" for item in files]

    def prepare(
        self,
        config: CaseConfig,
        destination: Path,
        *,
        generate_mesh: bool = False,
        reuse_mesh: bool = False,
        force: bool = False,
    ) -> Path:
        """Atomically produce a case directory, leaving no partial case on failure."""
        self._validate_options(generate_mesh, reuse_mesh)
        self._validate_base_case(config, reuse_mesh=reuse_mesh)
        destination = destination.resolve()
        if destination == self.base_case_directory:
            raise PreparationError(
                "Refusing to replace the base case; choose a separate output directory"
            )
        if destination.exists():
            is_nonempty_directory = destination.is_dir() and any(destination.iterdir())
            if (not destination.is_dir() or is_nonempty_directory) and not force:
                raise PreparationError(
                    f"Output '{destination}' already exists; pass --force to replace it"
                )

        destination.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(
            tempfile.mkdtemp(prefix=f".{destination.name}.staging-", dir=destination.parent)
        )
        try:
            self._materialize(
                config,
                staging,
                generate_mesh=generate_mesh,
                reuse_mesh=reuse_mesh,
            )
            if destination.exists():
                if destination.is_dir() and not destination.is_symlink():
                    shutil.rmtree(destination)
                else:
                    destination.unlink()
            staging.replace(destination)
        except Exception:
            shutil.rmtree(staging, ignore_errors=True)
            raise
        return destination

    def _materialize(
        self,
        config: CaseConfig,
        destination: Path,
        *,
        generate_mesh: bool,
        reuse_mesh: bool,
    ) -> None:
        self._copy_static_case_files(destination)
        self._write_text(
            destination / "constant" / "couplingProperties", render_coupling_properties(config)
        )

        indices = elmer_boundary_indices(len(config.electrodes.pairs))
        self._write_text(destination / "case.sif", self.elmer_renderer.render(config, indices))

        zero_directory = destination / "0"
        zero_directory.mkdir(parents=True, exist_ok=True)
        for name, content in self.openfoam_renderer.render(config).items():
            self._write_text(zero_directory / name, content)

        with (destination / "resolved-config.yaml").open(
            "w", encoding="utf-8", newline="\n"
        ) as stream:
            yaml.safe_dump(config.to_dict(), stream, sort_keys=False)

        if generate_mesh:
            self.mesh_generator.generate(config, destination / "channel.msh")
        elif reuse_mesh:
            source_mesh = self.base_case_directory / "channel.msh"
            if not source_mesh.is_file():
                raise PreparationError(f"Reusable base mesh not found: {source_mesh}")
            self._validate_reusable_mesh(config)
            shutil.copy2(source_mesh, destination / "channel.msh")

        self._write_manifest(config, destination, generate_mesh, reuse_mesh)

    def _copy_static_case_files(self, destination: Path) -> None:
        for directory_name in ("constant", "system"):
            source = self.base_case_directory / directory_name
            if not source.is_dir():
                raise PreparationError(f"Base case directory is missing '{source}'")
            shutil.copytree(
                source,
                destination / directory_name,
                dirs_exist_ok=True,
                ignore=(
                    shutil.ignore_patterns("polyMesh")
                    if directory_name == "constant"
                    else None
                ),
            )

    def _validate_base_case(self, config: CaseConfig, *, reuse_mesh: bool) -> None:
        for directory_name in ("constant", "system"):
            source = self.base_case_directory / directory_name
            if not source.is_dir():
                raise PreparationError(f"Base case directory is missing '{source}'")
        if reuse_mesh:
            source_mesh = self.base_case_directory / "channel.msh"
            if not source_mesh.is_file():
                raise PreparationError(f"Reusable base mesh not found: {source_mesh}")
            self._validate_reusable_mesh(config)

    def _write_manifest(
        self,
        config: CaseConfig,
        destination: Path,
        generated_mesh: bool,
        reused_mesh: bool,
    ) -> None:
        files = sorted(
            str(item.relative_to(destination))
            for item in destination.rglob("*")
            if item.is_file()
        )
        mesh_status = "omitted"
        if generated_mesh:
            mesh_status = "generated"
        elif reused_mesh:
            mesh_status = "reused"
        manifest: Dict[str, object] = {
            "schema_version": 1,
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "source_commit": self._source_commit(),
            "source_dirty": self._source_is_dirty(),
            "mesh": mesh_status,
            "files": files,
            "config": config.to_dict(),
        }
        self._write_text(destination / "manifest.json", json.dumps(manifest, indent=2) + "\n")

    def _source_commit(self) -> Optional[str]:
        result = subprocess.run(
            ("git", "rev-parse", "HEAD"),
            cwd=self.base_case_directory,
            check=False,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip() if result.returncode == 0 else None

    def _source_is_dirty(self) -> Optional[bool]:
        result = subprocess.run(
            ("git", "status", "--porcelain", "--untracked-files=normal"),
            cwd=self.base_case_directory,
            check=False,
            capture_output=True,
            text=True,
        )
        return bool(result.stdout.strip()) if result.returncode == 0 else None

    @staticmethod
    def _write_text(destination: Path, content: str) -> None:
        with destination.open("w", encoding="utf-8", newline="\n") as stream:
            stream.write(content)

    def _validate_reusable_mesh(self, config: CaseConfig) -> None:
        base_config_path = self.base_case_directory / "electrodes.yaml"
        try:
            base_config = load_case_config(base_config_path)
        except ConfigError as exc:
            raise PreparationError(
                f"Cannot validate the reusable mesh against '{base_config_path}': {exc}"
            ) from exc

        def signature(case: CaseConfig) -> object:
            return (
                case.channel,
                case.mesh,
                case.electrodes.length,
                tuple(pair.x_center for pair in case.electrodes.pairs),
            )

        if signature(config) != signature(base_config):
            raise PreparationError(
                "--reuse-mesh is only valid when channel dimensions, mesh controls, "
                "electrode count, centers, and shared length match the base "
                "electrodes.yaml; use --mesh"
            )

    @staticmethod
    def _validate_options(generate_mesh: bool, reuse_mesh: bool) -> None:
        if generate_mesh and reuse_mesh:
            raise PreparationError("Choose either --mesh or --reuse-mesh, not both")
