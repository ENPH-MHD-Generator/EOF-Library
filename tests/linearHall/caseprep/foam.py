"""Generic, deterministic rendering for OpenFOAM dictionary syntax."""

from __future__ import annotations

from pathlib import PurePosixPath
from typing import List, Mapping, Sequence, Union

from .models import OpenFoamFileConfig


FoamScalar = Union[None, bool, int, float, str]
FoamValue = Union[FoamScalar, Sequence["FoamValue"], Mapping[str, "FoamValue"]]


OPENFOAM_FILE_HEADER = """\
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  dev                                   |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {foam_class};
    location    "{location}";
    object      {object_name};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""

OPENFOAM_FILE_FOOTER = (
    "\n// ************************************************************************* //\n"
)


def _render_atom(value: FoamValue) -> str:
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return repr(value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, str):
        return value
    if isinstance(value, Mapping):
        return "{ " + " ".join(line.strip() for line in _render_entries(value)) + " }"
    if isinstance(value, Sequence):
        return "(" + " ".join(_render_atom(item) for item in value) + ")"
    raise TypeError(f"OpenFOAM dictionary value cannot be rendered inline: {value!r}")


def _render_entries(
    entries: Mapping[str, FoamValue], indentation: int = 0
) -> List[str]:
    prefix = " " * indentation
    lines: List[str] = []
    for key, value in entries.items():
        if isinstance(value, Mapping):
            lines.extend((f"{prefix}{key}", f"{prefix}{{"))
            lines.extend(_render_entries(value, indentation + 4))
            lines.append(f"{prefix}}}")
        elif value is None:
            lines.append(f"{prefix}{key};")
        else:
            lines.append(f"{prefix}{key:<20} {_render_atom(value)};")
    return lines


def render_openfoam_header(
    foam_class: str, location: str, object_name: str
) -> str:
    """Render the common FoamFile banner and metadata dictionary."""
    return OPENFOAM_FILE_HEADER.format(
        foam_class=foam_class,
        location=location,
        object_name=object_name,
    )


def render_openfoam_file(path: str, config: OpenFoamFileConfig) -> str:
    """Render one validated file using its path for location and object metadata."""
    parsed_path = PurePosixPath(path)
    header = render_openfoam_header(
        config.foam_class, parsed_path.parent.as_posix(), parsed_path.name
    )
    body = "\n".join(_render_entries(config.entries))
    return header + body + OPENFOAM_FILE_FOOTER
