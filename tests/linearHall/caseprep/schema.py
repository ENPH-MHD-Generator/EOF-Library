"""Generate the editor-facing JSON Schema from the Pydantic source models."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

from .models import CaseConfigV1


JSON_SCHEMA_DIALECT = "https://json-schema.org/draft/2020-12/schema"


def experiment_schema() -> Dict[str, Any]:
    """Return the JSON Schema for the current public experiment format."""
    generated = CaseConfigV1.model_json_schema(mode="validation")
    return {"$schema": JSON_SCHEMA_DIALECT, **generated}


def render_experiment_schema() -> str:
    return json.dumps(experiment_schema(), indent=2) + "\n"


def main(arguments: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate or verify the MHD experiment JSON Schema."
    )
    parser.add_argument("output", type=Path, help="Destination JSON Schema file.")
    parser.add_argument(
        "--check",
        action="store_true",
        help="Fail if the destination does not match the generated schema.",
    )
    args = parser.parse_args(arguments)

    rendered = render_experiment_schema()
    if args.check:
        try:
            current = args.output.read_text(encoding="utf-8")
        except OSError as exc:
            parser.exit(1, f"error: cannot read '{args.output}': {exc}\n")
        if current != rendered:
            parser.exit(1, f"error: generated schema is stale: {args.output}\n")
        return 0

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

