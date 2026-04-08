#!/usr/bin/env python3
"""Run parametric sweeps over electrode configurations.

Reads a study YAML file that references a base electrodes.yaml and defines
one or more sweeps.  For each case in each sweep, the script:

  1. Deep-copies the base config
  2. Applies the parameter override
  3. Creates an isolated case directory
  4. Runs configure.py to generate case.sif + OpenFOAM 0/ files
  5. Optionally regenerates the mesh and/or launches the simulation

Study YAML format
-----------------

    base_config: electrodes.yaml
    output_dir: "studies/{sweep}/{index:03d}"

    sweeps:
      # Simple single-parameter sweep
      - name: "pair_count"
        parameter: "electrodes.count"
        values: [2, 4, 6, 8]

      # Explicit list of case overrides (deep-merged into base)
      - name: "custom_layouts"
        cases:
          - electrodes:
              pairs:
                - { x_center: 0.05, length: 0.01, resistance: 1.0 }
                - { x_center: 0.15, length: 0.01, resistance: 1.0 }
          - electrodes:
              pairs:
                - { x_center: 0.04, length: 0.01, resistance: 1.0 }
                - { x_center: 0.08, length: 0.01, resistance: 1.0 }
                - { x_center: 0.12, length: 0.01, resistance: 1.0 }
                - { x_center: 0.16, length: 0.01, resistance: 1.0 }

Usage:
    python run_study.py study.yaml
    python run_study.py study.yaml --mesh           # regenerate mesh per case
    python run_study.py study.yaml --dry-run        # preview only
"""

import argparse
import json
import shutil
import subprocess
import sys
from copy import deepcopy
from pathlib import Path

import yaml


def _deep_merge(base, override):
    """Recursively merge override into base, returning a new dict."""
    result = deepcopy(base)
    for key, val in override.items():
        if isinstance(val, dict) and isinstance(result.get(key), dict):
            result[key] = _deep_merge(result[key], val)
        else:
            result[key] = deepcopy(val)
    return result


def _set_nested(d, dotted_path, value):
    """Set a value in a nested dict using a dotted key path.

    When the path contains 'electrodes.count', we also remove
    'electrodes.pairs' (and vice-versa) so that configure.py
    picks the intended mode unambiguously.
    """
    keys = dotted_path.split(".")
    target = d
    for key in keys[:-1]:
        target = target.setdefault(key, {})
    target[keys[-1]] = value

    # Resolve ambiguity between explicit pairs and count mode
    if keys[:1] == ["electrodes"]:
        elec = d.get("electrodes", {})
        if keys[-1] == "count":
            elec.pop("pairs", None)
        elif keys[-1] == "pairs":
            elec.pop("count", None)


def _expand_sweep(sweep, base_config):
    """Yield (index, case_config) for each case in the sweep."""
    if "cases" in sweep:
        for i, override in enumerate(sweep["cases"]):
            yield i, _deep_merge(base_config, override)
    elif "parameter" in sweep and "values" in sweep:
        for i, val in enumerate(sweep["values"]):
            cfg = deepcopy(base_config)
            _set_nested(cfg, sweep["parameter"], val)
            yield i, cfg
    else:
        raise ValueError(
            f"Sweep '{sweep.get('name', '?')}' must have 'cases' or "
            "'parameter'+'values'"
        )


def run_study(study_path, do_mesh=False, dry_run=False):
    study_dir = Path(study_path).parent
    with open(study_path) as f:
        study = yaml.safe_load(f)

    base_path = study_dir / study.get("base_config", "electrodes.yaml")
    with open(base_path) as f:
        base_config = yaml.safe_load(f)

    output_template = study.get("output_dir", "studies/{sweep}/{index:03d}")
    sweeps = study.get("sweeps", [])
    if not sweeps:
        print("No sweeps defined — nothing to do.")
        return

    manifest = []

    for sweep in sweeps:
        sweep_name = sweep.get("name", "sweep")
        print(f"\n{'='*60}")
        print(f"Sweep: {sweep_name}")
        print(f"{'='*60}")

        for index, case_config in _expand_sweep(sweep, base_config):
            case_dir = study_dir / output_template.format(
                sweep=sweep_name, index=index,
            )
            label = f"  [{sweep_name}/{index}]"

            if dry_run:
                print(f"{label} -> {case_dir}  (dry-run)")
                continue

            case_dir.mkdir(parents=True, exist_ok=True)

            # Write the resolved config for reproducibility
            case_yaml = case_dir / "electrodes.yaml"
            with open(case_yaml, "w") as f:
                yaml.dump(case_config, f, default_flow_style=False, sort_keys=False)

            # Copy static OpenFOAM directories (system, constant) if present
            for subdir in ("system", "constant"):
                src = study_dir / subdir
                dst = case_dir / subdir
                if src.is_dir() and not dst.exists():
                    shutil.copytree(src, dst)

            # Run configure.py in the case directory
            cmd = [
                sys.executable,
                str(study_dir / "configure.py"),
                "--config", str(case_yaml),
                "--output-dir", str(case_dir),
            ]
            if do_mesh:
                cmd.append("--mesh")

            print(f"{label} generating in {case_dir}")
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"{label} FAILED:")
                print(result.stderr)
                continue
            print(result.stdout.rstrip())

            manifest.append({
                "sweep": sweep_name,
                "index": index,
                "dir": str(case_dir),
                "config": case_config,
            })

    # Write manifest for post-processing
    if manifest and not dry_run:
        manifest_path = study_dir / "study_manifest.json"
        with open(manifest_path, "w") as f:
            json.dump(manifest, f, indent=2, default=str)
        print(f"\nManifest written to {manifest_path}")

    print("\nStudy complete.")


def main():
    parser = argparse.ArgumentParser(
        description="Run parametric electrode-placement sweeps.",
    )
    parser.add_argument(
        "study",
        help="Path to the study YAML file.",
    )
    parser.add_argument(
        "--mesh",
        action="store_true",
        help="Regenerate the Gmsh mesh for each case.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview cases without generating files.",
    )
    args = parser.parse_args()
    run_study(args.study, args.mesh, args.dry_run)


if __name__ == "__main__":
    main()
