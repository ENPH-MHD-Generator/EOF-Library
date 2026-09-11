"""Container command-line interface for preparing and running MHD cases."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

from .config import ConfigError
from .preparer import PreparationError
from .runtime import MhdRuntime


BASE_CASE_DIRECTORY = Path("/home/openfoam/EOF-Library/tests/linearHall")
EXPERIMENTS_DIRECTORY = Path("/experiments")
RUNS_DIRECTORY = Path("/runs")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mhd",
        description="Prepare and run configuration-driven MHD simulations.",
    )
    commands = parser.add_subparsers(dest="command", required=True)

    prepare = commands.add_parser(
        "prepare", help="Compile and convert one YAML case beneath /runs."
    )
    prepare.add_argument(
        "config", nargs="?", type=Path, help="YAML experiment configuration."
    )
    prepare.add_argument(
        "-l",
        "--list",
        action="store_true",
        help="List experiment YAML files available beneath /experiments.",
    )
    prepare.add_argument(
        "--name",
        help="Case directory name beneath /runs (default: configuration filename).",
    )
    prepare.add_argument(
        "--mesh",
        action="store_true",
        help="Generate a fresh mesh; the compatible base mesh is reused by default.",
    )
    prepare.add_argument(
        "--ranks", type=int, default=2, help="MPI ranks per solver (default: 2)."
    )
    prepare.add_argument(
        "--force", action="store_true", help="Replace an existing case with the same name."
    )
    prepare.add_argument(
        "--dry-run", action="store_true", help="Validate and show the preparation plan."
    )

    run = commands.add_parser("run", help="Execute one prepared case beneath /runs.")
    run.add_argument("name", nargs="?", help="Prepared case name beneath /runs.")
    run.add_argument(
        "-l",
        "--list",
        action="store_true",
        help="List valid prepared cases available beneath /runs.",
    )
    run.add_argument(
        "--no-postprocess",
        action="store_true",
        help="Skip reconstructPar and foamToVTK after a successful simulation.",
    )
    run.add_argument(
        "--dry-run", action="store_true", help="Validate and show the execution plan."
    )
    return parser


def main(
    arguments: Optional[Sequence[str]] = None,
    *,
    base_case_directory: Path = BASE_CASE_DIRECTORY,
    experiments_directory: Path = EXPERIMENTS_DIRECTORY,
    runs_directory: Path = RUNS_DIRECTORY,
) -> int:
    parser = build_parser()
    args = parser.parse_args(arguments)
    runtime = MhdRuntime(
        base_case_directory,
        runs_directory,
        experiments_directory=experiments_directory,
    )

    try:
        if args.command == "prepare":
            if args.list:
                if args.config is not None:
                    parser.error("prepare accepts either CONFIG or --list, not both")
                experiments = runtime.available_experiments()
                if not experiments:
                    print(
                        "No experiment YAML files found in "
                        f"{runtime.experiments_directory}"
                    )
                    return 0
                print("Available experiments:")
                for experiment in experiments:
                    relative = experiment.relative_to(runtime.experiments_directory)
                    print(f"  {relative}  (mhd prepare {experiment})")
                return 0
            if args.config is None:
                parser.error("prepare requires CONFIG or --list")
            destination = runtime.prepare(
                args.config,
                name=args.name,
                generate_mesh=args.mesh,
                ranks=args.ranks,
                force=args.force,
                dry_run=args.dry_run,
            )
            action = "Validated" if args.dry_run else "Prepared"
            print(f"{action} case: {destination}")
            return 0

        if args.list:
            if args.name is not None:
                parser.error("run accepts either NAME or --list, not both")
            cases = runtime.available_cases()
            if not cases:
                print(f"No runnable prepared cases found in {runtime.runs_directory}")
                return 0
            print("Available prepared cases:")
            for case_name in cases:
                print(f"  {case_name}  (mhd run {case_name})")
            return 0
        if args.name is None:
            parser.error("run requires NAME or --list")
        destination = runtime.run(
            args.name,
            postprocess=not args.no_postprocess,
            dry_run=args.dry_run,
        )
        action = "Validated" if args.dry_run else "Completed"
        print(f"{action} case: {destination}")
        return 0
    except (ConfigError, PreparationError, RuntimeError) as exc:
        parser.exit(2, f"error: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
