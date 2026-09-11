"""Shared boundary naming and tag conventions.

Gmsh, Elmer, and OpenFOAM must all consume this module. Keeping the mapping
here prevents a geometry implementation from silently changing the physical
meaning of a generated boundary condition.
"""

from typing import Dict, List


INLET = "InletX"
OUTLET = "OutletX"
INSULATOR = "InsulatorSurface"


def electrode_patch_names(pair_count: int) -> List[str]:
    """Return electrode patch names in deterministic pair order."""
    names: List[str] = []
    for index in range(1, pair_count + 1):
        names.extend((f"CathodeSurface_{index}", f"AnodeSurface_{index}"))
    return names


def boundary_tags(pair_count: int) -> Dict[str, int]:
    """Return the canonical Gmsh physical tags for every boundary."""
    tags = {INLET: 20, OUTLET: 21, INSULATOR: 30}
    for index in range(1, pair_count + 1):
        tags[f"CathodeSurface_{index}"] = 38 + 2 * index
        tags[f"AnodeSurface_{index}"] = 39 + 2 * index
    return tags


def elmer_boundary_indices(pair_count: int) -> Dict[str, int]:
    """Map boundary names to ElmerGrid's monotonically ordered indices."""
    tags = boundary_tags(pair_count)
    return {
        name: index
        for index, name in enumerate(sorted(tags, key=tags.get), start=1)
    }
