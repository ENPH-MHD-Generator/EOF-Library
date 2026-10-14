"""Tests for versioned, data-driven OpenFOAM case rendering."""

from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from pydantic import ValidationError

from caseprep.loader import load_case_config
from caseprep.models import CaseConfigV1, CaseConfigV2, OpenFoamConfig
from caseprep.preparer import CasePreparer
from caseprep.renderers import OpenFoamCaseRenderer


REPOSITORY = Path(__file__).resolve().parents[2]


class OpenFoamRenderingTest(TestCase):
    def test_current_experiment_loads_as_v2_and_renders_solver_files(self) -> None:
        config = load_case_config(REPOSITORY / "experiments" / "linear-hall.yaml")

        self.assertIsInstance(config, CaseConfigV2)
        rendered = OpenFoamCaseRenderer().render(config)
        self.assertIn(
            "rho                  0.325;",
            rendered["constant/transportProperties"],
        )
        self.assertIn("adjustTimeStep       yes;", rendered["system/controlDict"])
        self.assertIn("div(phi,U)", rendered["system/fvSchemes"])
        self.assertIn("$pcorr;", rendered["system/fvSolution"])

    def test_generic_renderer_supports_metadata_and_recursive_values(self) -> None:
        openfoam = OpenFoamConfig.model_validate(
            {
                "files": {
                    "constant/customState": {
                        "class": "uniformDimensionedVectorField",
                        "entries": {
                            "enabled": True,
                            "dimensions": "[0 1 -2 0 0 0 0]",
                            "value": [1, 2.5, 3],
                            "nested": {"$base": None, "mode": "cellPoint"},
                        },
                    }
                }
            }
        )
        config = load_case_config(REPOSITORY / "experiments" / "linear-hall.yaml")
        config = config.model_copy(update={"openfoam": openfoam})

        rendered = OpenFoamCaseRenderer().render(config)["constant/customState"]
        self.assertIn("class       uniformDimensionedVectorField;", rendered)
        self.assertIn('location    "constant";', rendered)
        self.assertIn("enabled              yes;", rendered)
        self.assertIn("value                (1 2.5 3);", rendered)
        self.assertIn("    $base;", rendered)

    def test_file_paths_cannot_escape_the_case(self) -> None:
        with self.assertRaisesRegex(ValidationError, "directly beneath"):
            OpenFoamConfig.model_validate(
                {"files": {"../system/controlDict": {"entries": {}}}}
            )

    def test_v1_remains_loadable_and_only_renders_initial_fields(self) -> None:
        config = load_case_config(REPOSITORY / "tests/linearHall/electrodes.yaml")

        self.assertIsInstance(config, CaseConfigV1)
        self.assertEqual(
            set(OpenFoamCaseRenderer().render(config)),
            {
                "0/U",
                "0/T",
                "0/p_rgh",
                "0/B",
                "0/Potential",
                "0/J_dens",
                "0/electric_field",
                "0/JH",
            },
        )

    def test_preparer_writes_configured_files_over_copied_baseline(self) -> None:
        config = load_case_config(REPOSITORY / "experiments" / "linear-hall.yaml")
        with TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            base = root / "base"
            (base / "constant").mkdir(parents=True)
            (base / "system").mkdir()
            (base / "constant/transportProperties").write_text(
                "stale\n", encoding="utf-8"
            )
            destination = root / "prepared"

            CasePreparer(base).prepare(config, destination)

            transport = (destination / "constant/transportProperties").read_text(
                encoding="utf-8"
            )
            self.assertNotIn("stale", transport)
            self.assertIn("rho                  0.325;", transport)
            self.assertTrue((destination / "system/fvSolution").is_file())

            resolved = load_case_config(destination / "resolved-config.yaml")
            self.assertEqual(resolved.to_dict(), config.to_dict())
