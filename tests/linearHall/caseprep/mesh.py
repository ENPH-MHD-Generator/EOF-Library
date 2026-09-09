"""Mesh generator adapters.

The Gmsh dependency is intentionally imported only when mesh generation is
requested. Configuration validation and text-file generation do not require it.
"""

from pathlib import Path

from .config import CaseConfig


class ProceduralMeshGenerator:
    """Generate the supported rectangular linear Hall channel mesh."""

    def generate(self, config: CaseConfig, destination: Path) -> None:
        try:
            from generate_mesh import generate
        except ModuleNotFoundError as exc:
            if exc.name == "gmsh":
                raise RuntimeError(
                    "Mesh generation requires the optional 'gmsh' Python package"
                ) from exc
            raise

        channel = config.channel
        mesh = config.mesh
        generate(
            out_msh=str(destination),
            mesh_size_min=mesh.size_min,
            mesh_size_max=mesh.size_max,
            mesh_size_factor=mesh.size_factor,
            channel_config={
                "num_pairs": len(config.electrodes.pairs),
                "channel_length": channel.length,
                "channel_height": channel.height,
                "channel_width": channel.width,
                "electrode_length": config.electrodes.length,
                "wall_thickness": channel.wall_thickness,
                "electrode_centers": [
                    pair.x_center for pair in config.electrodes.pairs
                ],
            },
        )
