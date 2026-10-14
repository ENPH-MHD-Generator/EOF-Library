"""Check that Saha physics and both mesh generators share the canonical models."""

import copy
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest import TestCase
from unittest.mock import Mock, patch

import yaml
from pydantic import ValidationError

from caseprep.loader import load_case_config
from caseprep.mesh import ProceduralMeshGenerator
from caseprep.models import CaseConfig
from caseprep.renderers import ElmerCaseRenderer, OpenFoamCaseRenderer, render_coupling_properties
from caseprep.boundaries import elmer_boundary_indices

ROOT = Path(__file__).resolve().parents[2]


class SahaMergeTest(TestCase):
    def setUp(self):
        self.raw = yaml.safe_load((ROOT / "experiments/linear-hall.yaml").read_text())

    def test_one_tetrahedral_size_policy_reaches_the_generator(self):
        self.raw["electrodes"] = {"count": 2, "length": 0.01}
        self.raw["mesh"]["target_element_size"] = 0.0123
        config = CaseConfig.model_validate(self.raw)
        tetra, structured = Mock(), Mock()
        with patch.dict(sys.modules, generate_mesh=SimpleNamespace(generate=tetra, generate_structured=structured)):
            ProceduralMeshGenerator().generate(config, Path("mesh.msh"))
        self.assertEqual(tetra.call_args.kwargs["target_element_size"], 0.0123)
        self.assertEqual(tetra.call_args.kwargs["channel_config"]["num_pairs"], 2)
        self.assertEqual(len(tetra.call_args.kwargs["channel_config"]["electrode_centers"]), 2)
        self.assertFalse(any(key.startswith("mesh_size") for key in tetra.call_args.kwargs))
        structured.assert_not_called()

    def test_structured_grading_reaches_its_generator_and_roundtrips(self):
        self.raw["mesh"] = {"type": "structured", "electrode_wall_cell_size": None, "side_wall_cell_size": None, "streamwise_cell_size": None}
        config = CaseConfig.model_validate(self.raw)
        self.assertEqual(CaseConfig.model_validate(config.to_dict()).to_dict(), config.to_dict())
        tetra, structured = Mock(), Mock()
        with patch.dict(sys.modules, generate_mesh=SimpleNamespace(generate=tetra, generate_structured=structured)):
            ProceduralMeshGenerator().generate(config, Path("mesh.msh"))
        self.assertIsNone(structured.call_args.kwargs["electrode_wall_cell_size"])
        self.assertIsNone(structured.call_args.kwargs["streamwise_cell_size"])
        tetra.assert_not_called()
        sif = ElmerCaseRenderer().render(config, elmer_boundary_indices(4))
        self.assertIn("Linear System Direct Method = MUMPS", sif)

    def test_mesh_types_reject_each_others_controls(self):
        for mesh in (
            {"type": "tetrahedral", "target_element_size": 0.005, "wall_cell_size": 0.001},
            {"type": "structured", "target_element_size": 0.005},
            {"type": "tetrahedral", "target_element_size": 0.005, "size_factor": 0.35},
        ):
            with self.subTest(mesh=mesh), self.assertRaises(ValidationError):
                CaseConfig.model_validate(dict(self.raw, mesh=mesh))

    def test_structured_cross_field_constraints(self):
        for settings in (
            {"growth_rate": 1}, {"wall_cell_size": 0.01},
            {"side_wall_cell_size": 0.01}, {"electrode_wall_cell_size": 0.01},
            {"electrode_edge_cell_size": 0.02},
        ):
            with self.subTest(settings=settings), self.assertRaises(ValidationError):
                CaseConfig.model_validate(dict(self.raw, mesh={"type": "structured", **settings}))

    def test_plasma_bounds_and_strict_switches(self):
        for changes in (
            {"seed_mole_fraction": 1}, {"sigma_min": 1e7},
            {"electron_temperature_relaxation": 1.1}, {"two_temperature": "false"},
            {"electron_energy_transport": 1}, {"potassium_elastic_scale": float("nan")},
            {"electron_transport_model": "unknown"}, {"seed_cross_section": 1e-19},
        ):
            with self.subTest(changes=changes), self.assertRaises(ValidationError):
                CaseConfig.model_validate(dict(self.raw, plasma={**self.raw["plasma"], **changes}))

    def test_coupling_requires_nonnegative_integer_step_count(self):
        for value in (-1, 1.5, True):
            with self.subTest(value=value), self.assertRaises(ValidationError):
                CaseConfig.model_validate(dict(self.raw, coupling={"max_steps_between_updates": value}))

    def test_rejects_incompressible_and_duplicate_coupling_inputs(self):
        for path in ("constant/transportProperties", "0/p_rgh", "constant/couplingProperties"):
            raw = copy.deepcopy(self.raw)
            raw["openfoam"]["files"][path] = {"entries": {}}
            with self.subTest(path=path), self.assertRaises(ValidationError):
                CaseConfig.model_validate(raw)

    def test_saha_rendering_preserves_physics_and_resolved_auto_pairs(self):
        self.raw["electrodes"] = {"count": 2, "length": 0.01}
        self.raw["physics"]["insulator_wall_temperature"] = None
        self.raw["physics"]["electrode_wall_temperature"] = 450
        config = CaseConfig.model_validate(self.raw)
        rendered = OpenFoamCaseRenderer().render(config)
        self.assertIn("0/p", rendered)
        self.assertIn("0/Te", rendered)
        self.assertNotIn("0/p_rgh", rendered)
        self.assertIn("uniform 450", rendered["0/T"])
        self.assertIn("CathodeSurface_2", rendered["0/T"])
        self.assertNotIn("CathodeSurface_3", rendered["0/T"])
        sif = ElmerCaseRenderer().render(config, elmer_boundary_indices(2), {"ionization_energy": 4.5, "weight_ratio": 0.75})
        self.assertIn("ElectronEnergySolver", sif)
        self.assertIn("4.5", sif)
        self.assertIn("Element Joule Heating", sif)
        self.assertIn("maxStepsBetweenUpdates", render_coupling_properties(config))
