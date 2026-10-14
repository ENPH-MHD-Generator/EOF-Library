"""Regression coverage for canonical configs and the prepare/run lifecycle."""

import json
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase
from unittest.mock import Mock, patch

import yaml
from pydantic import ValidationError

from caseprep.cli import main
from caseprep.foam import render_openfoam_file
from caseprep.loader import ConfigError, load_case_config
from caseprep.models import CaseConfig, OpenFoamConfig, OpenFoamFileConfig
from caseprep.preparer import CasePreparer, PreparationError
from caseprep.runtime import MhdRuntime, RuntimeCommandError
from caseprep.schema import render_experiment_schema
from caseprep.collisions import TABLE_FILE_NAME, METADATA_FILE_NAME
from types import SimpleNamespace


def fake_collision_tables(config, destination):
    (destination / TABLE_FILE_NAME).write_text("test collision data\n")
    (destination / METADATA_FILE_NAME).write_text("{}\n")
    return SimpleNamespace(saha={"ionization_energy": 4.341, "weight_ratio": 0.5})


REPOSITORY = Path(__file__).resolve().parents[2]
BASE = REPOSITORY / "tests/linearHall"
EXPERIMENT = REPOSITORY / "experiments/linear-hall.yaml"


class CanonicalConfigTest(TestCase):
    def setUp(self) -> None:
        self.raw = yaml.safe_load(EXPERIMENT.read_text())

    def test_the_only_schema_artifact_matches_the_model(self) -> None:
        schemas = list((REPOSITORY / "schemas").glob("*.schema.json"))
        self.assertEqual([path.name for path in schemas], ["mhd-experiment.schema.json"])
        self.assertEqual(schemas[0].read_text(), render_experiment_schema())

    def test_invalid_geometry_and_numbers_are_rejected(self) -> None:
        for value in (True, "0.2", 0, -1, float("inf"), float("nan")):
            with self.subTest(length=value):
                self.raw["channel"]["length"] = value
                with self.assertRaises(ValidationError):
                    CaseConfig.model_validate(self.raw)
        self.raw["channel"]["length"] = 0.2
        for centers in ((0, 0.1), (0.05, 0.051), (0.1, 0.2)):
            with self.subTest(centers=centers):
                self.raw["electrodes"]["pairs"] = [
                    {"x_center": center, "resistance": 1} for center in centers
                ]
                with self.assertRaises(ValidationError):
                    CaseConfig.model_validate(self.raw)

    def test_automatic_placement_roundtrips_to_explicit_canonical_config(self) -> None:
        self.raw["electrodes"] = {"length": 0.01, "count": 4, "resistance": 2}
        config = CaseConfig.model_validate(self.raw)
        self.assertEqual(len(config.electrode_pairs), 4)
        self.assertEqual(CaseConfig.model_validate(config.to_dict()).to_dict(), config.to_dict())
        self.assertTrue(all(pair.resistance == 2 for pair in config.electrode_pairs))

    def test_mesh_rejects_unknown_sizing_controls(self) -> None:
        self.raw["mesh"] = {"type": "tetrahedral", "size_min": 0.005}
        with self.assertRaisesRegex(ValidationError, "target_element_size"):
            CaseConfig.model_validate(self.raw)

    def test_invalid_yaml_and_duplicate_keys_fail_clearly(self) -> None:
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / "config.yaml"
            for document, message in (
                ("{}", "Field required"),
                ("[]", "YAML mapping"),
                ("mesh: {}\nmesh: {}", "duplicate key"),
                ("channel: {length: 1, length: 2}", "duplicate key"),
            ):
                with self.subTest(document=document):
                    path.write_text(document)
                    with self.assertRaisesRegex(ConfigError, message):
                        load_case_config(path)

    def test_noncanonical_paths_are_rejected(self) -> None:
        for path in ("system/../x", "/system/x", "system//x", "system/./x", "system/x/", "0/", "system/bad name", 'system/bad"name'):
            with self.subTest(path=path), self.assertRaises(ValidationError):
                OpenFoamConfig.model_validate({"files": {path: {"entries": {}}}})

    def test_recursive_lists_render_and_preserve_numeric_precision(self) -> None:
        config = OpenFoamFileConfig(entries={
            "values": [[1, 2], {"enabled": True, "$base": None}],
            "precise": 1.2345678901234567,
        })
        rendered = render_openfoam_file("system/custom", config)
        self.assertIn("((1 2) { enabled              yes; $base; })", rendered)
        self.assertIn("1.2345678901234567;", rendered)

    def test_nonfinite_values_and_null_list_items_fail_before_rendering(self) -> None:
        for value in (float("nan"), float("inf"), [None], {"nested": [float("inf")]}):
            with self.subTest(value=value), self.assertRaises(ValidationError):
                OpenFoamFileConfig(entries={"value": value})


class PipelineTest(TestCase):
    def setUp(self) -> None:
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.runtime = MhdRuntime(BASE, self.root / "runs", experiments_directory=REPOSITORY / "experiments")
        self.config = load_case_config(EXPERIMENT)
        self.collision_writer = self.enterContext(patch(
            "caseprep.preparer.write_collision_tables", side_effect=fake_collision_tables
        ))

    @staticmethod
    def converted_mesh(command, case_directory, **kwargs) -> None:
        """Stand in for native converters, retaining their actual output layout."""
        if command[0] == "gmshToFoam":
            mesh = case_directory / "constant/polyMesh"
            mesh.mkdir(parents=True)
            for filename in ("points", "faces", "owner", "neighbour", "boundary"):
                (mesh / filename).write_text("mesh\n")
        elif command[0] == "decomposePar":
            import shutil
            import re
            content = (case_directory / "system/decomposeParDict").read_text()
            ranks = int(re.search(r"numberOfSubdomains\s+(\d+)", content).group(1))
            for index in range(ranks):
                shutil.copytree(case_directory / "constant/polyMesh", case_directory / f"processor{index}/constant/polyMesh")
        elif command[0] == "ElmerGrid":
            mesh = case_directory / "meshElmer"
            mesh.mkdir(exist_ok=True)
            for filename in ("mesh.header", "mesh.nodes", "mesh.elements", "mesh.boundary"):
                (mesh / filename).write_text("mesh\n")
            if "-metis" in command:
                ranks = int(command[-1])
                if ranks == 1:
                    return
                partition = mesh / f"partitioning.{ranks}"
                partition.mkdir()
                for index in range(1, ranks + 1):
                    for suffix in ("header", "nodes", "elements", "boundary", "shared"):
                        (partition / f"part.{index}.{suffix}").write_text("mesh\n")

    def prepare(self) -> Path:
        with patch.object(self.runtime, "_run", side_effect=self.converted_mesh):
            return self.runtime.prepare(EXPERIMENT, name="case")

    def test_default_example_reuses_mesh_and_writes_canonical_marker(self) -> None:
        destination = self.prepare()
        marker = json.loads((destination / self.runtime.MARKER_NAME).read_text())
        self.assertEqual(marker["config"], self.config.to_dict())
        self.assertEqual(marker["ranks"], 2)
        self.assertEqual((destination / "channel.msh").read_bytes(), (BASE / "channel.msh").read_bytes())
        self.assertEqual(self.runtime.available_cases(), ["case"])

    def test_coarse_example_requires_new_mesh(self) -> None:
        config = load_case_config(REPOSITORY / "experiments/linear-hall-coarse.yaml")
        with self.assertRaisesRegex(PreparationError, "use --mesh"):
            CasePreparer(BASE).describe(config, self.root / "case", reuse_mesh=True)

    def test_failed_native_command_preserves_existing_case(self) -> None:
        destination = self.prepare()
        marker = (destination / self.runtime.MARKER_NAME).read_bytes()
        with patch.object(self.runtime, "_run", side_effect=RuntimeCommandError("conversion failed")):
            with self.assertRaisesRegex(RuntimeCommandError, "conversion failed"):
                self.runtime.prepare(EXPERIMENT, name="case", force=True)
        self.assertEqual((destination / self.runtime.MARKER_NAME).read_bytes(), marker)
        self.assertEqual(list(self.runtime.runs_directory.iterdir()), [destination])

    def test_run_uses_recorded_ranks_and_optional_postprocessing(self) -> None:
        self.prepare()
        with patch.object(self.runtime, "_run") as run, patch.object(self.runtime, "_postprocess") as postprocess:
            self.runtime.run("case")
        commands = [call.args[0] for call in run.call_args_list]
        self.assertEqual(commands[0][-9:], ["-n", "2", "mdhLinearHall", "-parallel", ":", "-n", "2", "ElmerSolver_mpi", "case.sif"])
        postprocess.assert_called_once()
        with patch.object(self.runtime, "_run") as run:
            self.runtime.run("case", postprocess=False)
            self.assertEqual(run.call_count, 1)

    def test_missing_partition_is_rejected_by_listing_and_run(self) -> None:
        destination = self.prepare()
        (destination / "meshElmer/partitioning.2/part.1.header").unlink()
        self.assertEqual(self.runtime.available_cases(), [])
        with self.assertRaisesRegex(PreparationError, "part.1.header"):
            self.runtime.run("case", dry_run=True)

    def test_single_rank_marker_is_rejected_before_launch(self) -> None:
        destination = self.prepare()
        marker_path = destination / self.runtime.MARKER_NAME
        marker = json.loads(marker_path.read_text())
        marker["ranks"] = 1
        marker_path.write_text(json.dumps(marker))
        self.assertEqual(self.runtime.available_cases(), [])
        with patch.object(self.runtime, "_run") as run:
            with self.assertRaisesRegex(PreparationError, "at least 2"):
                self.runtime.run("case")
            run.assert_not_called()

    def test_missing_openfoam_mesh_file_is_rejected(self) -> None:
        destination = self.prepare()
        (destination / "processor0/constant/polyMesh/points").unlink()
        self.assertEqual(self.runtime.available_cases(), [])
        with self.assertRaisesRegex(PreparationError, "polyMesh/points"):
            self.runtime.run("case", dry_run=True)

    def test_invalid_names_and_ranks_are_rejected(self) -> None:
        for name in ("../escape", "/absolute", ".", "a/b"):
            with self.subTest(name=name), self.assertRaises(PreparationError):
                self.runtime.prepare(EXPERIMENT, name=name, dry_run=True)
        for ranks in (True, 0, 1, -1, 1.5):
            with self.subTest(ranks=ranks), self.assertRaises(PreparationError):
                self.runtime.prepare(EXPERIMENT, ranks=ranks, dry_run=True)

    def test_cli_dry_run_validates_without_writing(self) -> None:
        with redirect_stdout(StringIO()):
            result = main(["prepare", str(EXPERIMENT), "--dry-run"], base_case_directory=BASE, runs_directory=self.runtime.runs_directory)
        self.assertEqual(result, 0)
        self.assertFalse(self.runtime.runs_directory.exists())

    def test_preparer_restores_previous_directory_if_promotion_fails(self) -> None:
        destination = self.root / "case"
        destination.mkdir()
        (destination / "keep").write_text("previous")
        original_replace = Path.replace

        def fail_promotion(path, target):
            if ".staging-" in path.name:
                raise OSError("promotion failed")
            return original_replace(path, target)

        with patch.object(Path, "replace", fail_promotion):
            with self.assertRaisesRegex(OSError, "promotion failed"):
                CasePreparer(BASE).prepare(self.config, destination, force=True)
        self.assertEqual((destination / "keep").read_text(), "previous")
        self.assertEqual(list(self.root.iterdir()), [destination])

    def test_parallel_vtk_export_uses_time_names_to_avoid_overwriting_results(self) -> None:
        self.root.mkdir(exist_ok=True)
        process = Mock(returncode=0)
        with patch("caseprep.runtime.subprocess.Popen", return_value=process) as popen, patch("caseprep.runtime.shutil.which", return_value="/usr/bin/foamToVTK"):
            self.runtime._parallel_over_times("foamToVTK", ["0", "1e-8", "2e-8"], 2, self.root, [])
        self.assertEqual(popen.call_count, 2)
        for call in popen.call_args_list:
            self.assertIn("-useTimeName", call.args[0])

    def test_invalid_prepared_config_is_rejected_before_launch(self) -> None:
        destination = self.prepare()
        marker_path = destination / self.runtime.MARKER_NAME
        marker = json.loads(marker_path.read_text())
        marker["config"]["mesh"] = {"type": "unknown"}
        marker_path.write_text(json.dumps(marker))
        self.assertEqual(self.runtime.available_cases(), [])
        with patch.object(self.runtime, "_run") as run:
            with self.assertRaisesRegex(PreparationError, "Invalid configuration"):
                self.runtime.run("case")
            run.assert_not_called()
