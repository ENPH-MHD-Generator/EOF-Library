"""Generate both supported meshes and verify their shared solver boundary contract."""

from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

import gmsh
import yaml

from caseprep.boundaries import boundary_tags
from caseprep.mesh import ProceduralMeshGenerator
from caseprep.models import CaseConfig

ROOT = Path(__file__).resolve().parents[2]


class MeshGenerationTest(TestCase):
    def test_tetrahedral_and_structured_meshes_share_named_boundaries(self):
        raw = yaml.safe_load((ROOT / "experiments/linear-hall.yaml").read_text())
        raw["channel"] = {"length": 0.04, "height": 0.02, "width": 0.02, "wall_thickness": 0.002}
        raw["electrodes"] = {"length": 0.005, "count": 1, "resistance": 1}
        meshes = (
            ({"type": "tetrahedral", "target_element_size": 0.01}, 4),
            ({"type": "structured", "cell_size": 0.01, "wall_cell_size": 0.002, "electrode_wall_cell_size": None, "side_wall_cell_size": None, "streamwise_cell_size": 0.01, "electrode_edge_cell_size": 0.002, "growth_rate": 1.5}, 5),
        )
        with TemporaryDirectory() as temporary:
            for mesh, expected_type in meshes:
                with self.subTest(mesh=mesh["type"]):
                    config = CaseConfig.model_validate(dict(raw, mesh=mesh))
                    path = Path(temporary) / (mesh["type"] + ".msh")
                    ProceduralMeshGenerator().generate(config, path)
                    gmsh.initialize()
                    try:
                        gmsh.option.setNumber("General.Terminal", 0)
                        gmsh.open(str(path))
                        names = {gmsh.model.getPhysicalName(dim, tag) for dim, tag in gmsh.model.getPhysicalGroups()}
                        self.assertEqual(names, set(boundary_tags(1)) | {"Plasma"})
                        types, elements, _ = gmsh.model.mesh.getElements(3)
                        self.assertEqual(set(types), {expected_type})
                        self.assertTrue(all(len(items) > 0 for items in elements))
                        for _, tag in gmsh.model.getPhysicalGroups(2):
                            self.assertTrue(len(gmsh.model.getEntitiesForPhysicalGroup(2, tag)) > 0)
                    finally:
                        gmsh.finalize()
