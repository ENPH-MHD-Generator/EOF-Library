"""Mesh generator adapters.

The Gmsh dependency is intentionally imported only when mesh generation is
requested. Configuration validation and text-file generation do not require it.
"""

from pathlib import Path

from .models import CaseConfig


class ProceduralMeshGenerator:
    """Generate the supported rectangular linear Hall channel mesh."""

    def generate(self, config: CaseConfig, destination: Path) -> None:
        try:
            from generate_mesh import generate, generate_structured
        except ModuleNotFoundError as exc:
            if exc.name == "gmsh":
                raise RuntimeError(
                    "Mesh generation requires the optional 'gmsh' Python package"
                ) from exc
            raise

        mesh = config.mesh
        if mesh.type == "structured":
            generate_structured(
                out_msh=str(destination),
                channel_config=self._channel_config(config),
                cell_size=mesh.cell_size,
                wall_cell_size=mesh.wall_cell_size,
                electrode_edge_cell_size=mesh.electrode_edge_cell_size,
                growth_rate=mesh.growth_rate,
                electrode_wall_cell_size=mesh.electrode_wall_cell_size,
                side_wall_cell_size=mesh.side_wall_cell_size,
                streamwise_cell_size=mesh.streamwise_cell_size,
            )
            return
        generate(
            out_msh=str(destination),
            target_element_size=mesh.target_element_size,
            channel_config=self._channel_config(config),
        )

    @staticmethod
    def _channel_config(config: CaseConfig) -> dict:
        channel = config.channel
        return {
            "num_pairs": len(config.electrode_pairs),
            "channel_length": channel.length,
            "channel_height": channel.height,
            "channel_width": channel.width,
            "electrode_length": config.electrodes.length,
            "wall_thickness": channel.wall_thickness,
            "electrode_centers": [pair.x_center for pair in config.electrode_pairs],
        }
