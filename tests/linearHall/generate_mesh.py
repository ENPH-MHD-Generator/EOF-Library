import argparse
import math
from pathlib import Path

import gmsh

from caseprep.boundaries import boundary_tags


# ─── Default geometry constants ────────────────────────────────────
DEFAULTS = {
    "num_pairs": 4,
    "channel_length": 0.200,    # m
    "channel_height": 0.050,    # m
    "channel_width": 0.050,     # m
    "electrode_length": 0.010,  # m
    "wall_thickness": 0.005,    # m
}

# ─── Material tags (3-D physical groups) ───────────────────────────
MATERIAL_TAGS = {
    "Plasma": 1,
    "Cathode": 2,
    "Anode": 3,
    "Insulator": 4,
}

MATERIAL_PRIORITY = {
    "Plasma": 0,
    "Insulator": 1,
    "Cathode": 2,
    "Anode": 3,
}


# ─── Public helpers used by the procedural mesh generator ──────────────

def build_boundary_tags(num_pairs):
    """Backward-compatible alias for the shared boundary convention."""
    return boundary_tags(num_pairs)


def compute_electrode_centers(channel_length, num_pairs, explicit_centers=None):
    """Return list of x-coordinates for electrode-pair centers."""
    if explicit_centers is not None:
        if len(explicit_centers) != num_pairs:
            raise ValueError(
                f"Expected {num_pairs} electrode centers, got {len(explicit_centers)}"
            )
        return list(explicit_centers)
    spacing = channel_length / (num_pairs + 1)
    return [(i + 1) * spacing for i in range(num_pairs)]


# ─── Internal helpers (unchanged logic) ────────────────────────────

def _entity_bbox(dim: int, tag: int):
    return gmsh.model.getBoundingBox(dim, tag)


def _collect_external_surfaces(volume_tags):
    surface_to_n_up = {}
    for vtag in volume_tags:
        for dim, stag in gmsh.model.getBoundary([(3, vtag)], oriented=False, recursive=False):
            if dim != 2:
                continue
            up, _ = gmsh.model.getAdjacencies(2, stag)
            surface_to_n_up[stag] = len(up)
    return sorted(stag for stag, n_up in surface_to_n_up.items() if n_up == 1)


def _verify_every_boundary_has_exactly_one_physical(volume_tags):
    ext_surfs = _collect_external_surfaces(volume_tags)
    missing = []
    multiply_tagged = []
    for stag in ext_surfs:
        phys = gmsh.model.getPhysicalGroupsForEntity(2, stag)
        if len(phys) == 0:
            missing.append(stag)
        elif len(phys) > 1:
            multiply_tagged.append((stag, phys))
    if missing or multiply_tagged:
        raise RuntimeError(
            "Boundary tagging validation failed. "
            f"Missing: {missing}; multiply tagged: {multiply_tagged}"
        )


def _verify_tetra_only():
    elem_types, _, _ = gmsh.model.mesh.getElements(3)
    bad_types = []
    for et in elem_types:
        name, dim, _, _, _, _ = gmsh.model.mesh.getElementProperties(et)
        if dim != 3 or "tetrahedron" not in name.lower():
            bad_types.append((et, name))
    if bad_types:
        raise RuntimeError(
            "3D mesh contains non-tetra elements: "
            + ", ".join(f"type={et}:{name}" for et, name in bad_types)
        )


def _verify_all_boundary_faces_mapped(volume_tags):
    ext_surfs = _collect_external_surfaces(volume_tags)
    for stag in ext_surfs:
        phys = gmsh.model.getPhysicalGroupsForEntity(2, stag)
        if len(phys) != 1:
            raise RuntimeError(
                f"Exterior surface {stag} has {len(phys)} physical groups; expected exactly 1."
            )


def _clear_existing_physical_groups():
    for dim, tag in gmsh.model.getPhysicalGroups():
        gmsh.model.removePhysicalGroups([(dim, tag)])


def _configure_mesh_sizing(mesh_size_min, mesh_size_max, mesh_size_factor):
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFactor", mesh_size_factor)
    if mesh_size_min is not None:
        gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size_min)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh_size_min)
    if mesh_size_max is not None:
        gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size_max)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size_max)


def _build_geometry_and_assign_materials(cfg, centers):
    """Build a rectangular Hall channel with evenly-spaced electrode pairs."""
    occ = gmsh.model.occ
    L = cfg["channel_length"]
    H = cfg["channel_height"]
    W = cfg["channel_width"]
    t = cfg["wall_thickness"]
    el = cfg["electrode_length"]

    plasma = occ.addBox(0, 0, 0, L, H, W)
    imported = [(3, plasma)]
    imported_materials = ["Plasma"]

    for box_args in [
        (0, -t, 0, L, t, W),    # bottom wall
        (0, H, 0, L, t, W),     # top wall
        (0, 0, -t, L, H, t),    # front wall (z = 0)
        (0, 0, W, L, H, t),     # back  wall (z = W)
    ]:
        tag = occ.addBox(*box_args)
        imported.append((3, tag))
        imported_materials.append("Insulator")

    for xc in centers:
        x0 = xc - el / 2
        cat = occ.addBox(x0, -t, 0, el, t, W)
        imported.append((3, cat))
        imported_materials.append("Cathode")
        ano = occ.addBox(x0, H, 0, el, t, W)
        imported.append((3, ano))
        imported_materials.append("Anode")

    _, out_map = occ.fragment(imported, [])
    occ.synchronize()

    material_vols = {name: [] for name in MATERIAL_TAGS}
    if len(out_map) != len(imported):
        raise RuntimeError(
            f"Unexpected OCC fragment map size: got {len(out_map)}, expected {len(imported)}."
        )

    for i, children in enumerate(out_map):
        material = imported_materials[i]
        for dim, tag in children:
            if dim == 3:
                material_vols[material].append(tag)

    all_model_vols = {tag for dim, tag in gmsh.model.getEntities(3) if dim == 3}
    vol_to_materials = {v: [] for v in all_model_vols}
    for material, vols in material_vols.items():
        for v in set(vols):
            if v in vol_to_materials:
                vol_to_materials[v].append(material)

    unassigned = sorted(v for v, mats in vol_to_materials.items() if len(mats) == 0)
    if unassigned:
        raise RuntimeError(
            f"Post-fragment material mapping is invalid. Unassigned volumes: {unassigned}"
        )

    resolved_vol_to_material = {}
    multiply_assigned = []
    for v, mats in vol_to_materials.items():
        unique_mats = sorted(set(mats), key=lambda m: MATERIAL_PRIORITY[m], reverse=True)
        if not unique_mats:
            continue
        if len(unique_mats) > 1:
            multiply_assigned.append((v, unique_mats))
        resolved_vol_to_material[v] = unique_mats[0]

    if multiply_assigned:
        print(
            "[gmsh] INFO: resolved multi-material fragment volumes by priority: "
            f"{multiply_assigned}"
        )

    material_vols = {name: [] for name in MATERIAL_TAGS}
    for v, material in resolved_vol_to_material.items():
        material_vols[material].append(v)

    for material_name, tag in MATERIAL_TAGS.items():
        vols = sorted(set(material_vols[material_name]))
        if not vols:
            raise RuntimeError(f"Material '{material_name}' has no assigned volume.")
        gmsh.model.addPhysicalGroup(3, vols, tag=tag)
        gmsh.model.setPhysicalName(3, tag, material_name)

    vol_to_material = {}
    for material, vols in material_vols.items():
        for v in vols:
            vol_to_material[v] = material

    new_vols = sorted(all_model_vols)
    if not new_vols:
        raise RuntimeError("No volumes present after OCC fragment/synchronize.")
    return sorted(new_vols), vol_to_material


def _collect_material_interface_surfaces(vol_to_material, num_pairs, centers):
    """Find interface surfaces between plasma and surrounding materials."""
    insulator_surfs = []
    cathode_surfs = []
    anode_surfs = []

    for _, stag in gmsh.model.getEntities(2):
        up, _ = gmsh.model.getAdjacencies(2, stag)
        if len(up) != 2:
            continue
        mats = {vol_to_material.get(up[0]), vol_to_material.get(up[1])}
        if mats == {"Plasma", "Insulator"}:
            insulator_surfs.append(stag)
        elif mats == {"Plasma", "Anode"}:
            anode_surfs.append(stag)
        elif mats == {"Plasma", "Cathode"}:
            cathode_surfs.append(stag)

    cathode_groups = {i: [] for i in range(num_pairs)}
    anode_groups = {i: [] for i in range(num_pairs)}

    for stag in cathode_surfs:
        xmin, _, _, xmax, _, _ = gmsh.model.getBoundingBox(2, stag)
        x_mid = (xmin + xmax) / 2
        pair_idx = min(range(num_pairs),
                       key=lambda j: abs(x_mid - centers[j]))
        cathode_groups[pair_idx].append(stag)

    for stag in anode_surfs:
        xmin, _, _, xmax, _, _ = gmsh.model.getBoundingBox(2, stag)
        x_mid = (xmin + xmax) / 2
        pair_idx = min(range(num_pairs),
                       key=lambda j: abs(x_mid - centers[j]))
        anode_groups[pair_idx].append(stag)

    return sorted(set(insulator_surfs)), cathode_groups, anode_groups


# ─── Public entry point ────────────────────────────────────────────

def generate(
    out_msh="channel.msh",
    mesh_size_min=None,
    mesh_size_max=None,
    mesh_size_factor=1.0,
    channel_config=None,
):
    """Generate the Gmsh mesh.

    Parameters
    ----------
    channel_config : dict, optional
        Keys: num_pairs, channel_length, channel_height, channel_width,
        electrode_length, wall_thickness, electrode_centers (list or None).
        Missing keys fall back to DEFAULTS.

    Returns
    -------
    boundary_tags : dict
        Mapping of boundary name -> physical group tag.
    """
    cfg = {**DEFAULTS, **(channel_config or {})}
    num_pairs = cfg["num_pairs"]
    centers = compute_electrode_centers(
        cfg["channel_length"], num_pairs, cfg.get("electrode_centers"),
    )
    boundary_tags = build_boundary_tags(num_pairs)

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1)
        _configure_mesh_sizing(mesh_size_min, mesh_size_max, mesh_size_factor)
        gmsh.model.add("linear_hall_channel")

        volume_tags, vol_to_material = _build_geometry_and_assign_materials(cfg, centers)

        ins_surfs, cathode_groups, anode_groups = _collect_material_interface_surfaces(
            vol_to_material, num_pairs, centers,
        )
        if not ins_surfs:
            raise RuntimeError("No Plasma-Insulator interface surfaces found.")

        interface_set = set(ins_surfs)
        for i in range(num_pairs):
            if not cathode_groups[i]:
                raise RuntimeError(f"No cathode surfaces found for electrode pair {i + 1}.")
            if not anode_groups[i]:
                raise RuntimeError(f"No anode surfaces found for electrode pair {i + 1}.")
            interface_set.update(cathode_groups[i])
            interface_set.update(anode_groups[i])

        solid_vols = [(3, v) for v in volume_tags if vol_to_material.get(v) != "Plasma"]
        if solid_vols:
            gmsh.model.occ.remove(solid_vols, recursive=False)
            gmsh.model.occ.synchronize()

        _clear_existing_physical_groups()
        plasma_vols = [tag for dim, tag in gmsh.model.getEntities(3) if dim == 3]
        if not plasma_vols:
            raise RuntimeError("No plasma volume remains after removing solids.")
        gmsh.model.addPhysicalGroup(3, sorted(plasma_vols), tag=MATERIAL_TAGS["Plasma"])
        gmsh.model.setPhysicalName(3, MATERIAL_TAGS["Plasma"], "Plasma")

        gxmin, gymin, gzmin, gxmax, gymax, gzmax = gmsh.model.getBoundingBox(-1, -1)
        dom = max(gxmax - gxmin, gymax - gymin, gzmax - gzmin)
        tol = max(1e-9, 1e-6 * dom)

        groups = {name: [] for name in boundary_tags}
        exterior_surfs = _collect_external_surfaces(plasma_vols)
        ins_set = set(ins_surfs)

        for stag in exterior_surfs:
            xmin, _, _, xmax, _, _ = _entity_bbox(2, stag)
            if abs(xmin - xmax) <= tol and abs(xmin - gxmin) <= tol:
                groups["InletX"].append(stag)
                continue
            if abs(xmin - xmax) <= tol and abs(xmax - gxmax) <= tol:
                groups["OutletX"].append(stag)
                continue

            if stag in ins_set:
                groups["InsulatorSurface"].append(stag)
                continue

            matched = False
            for i in range(num_pairs):
                if stag in cathode_groups[i]:
                    groups[f"CathodeSurface_{i + 1}"].append(stag)
                    matched = True
                    break
                if stag in anode_groups[i]:
                    groups[f"AnodeSurface_{i + 1}"].append(stag)
                    matched = True
                    break
            if not matched:
                raise RuntimeError(
                    f"Plasma exterior surface not classified: {stag}"
                )

        for bname, ptag in boundary_tags.items():
            stags = sorted(set(groups[bname]))
            if not stags:
                raise RuntimeError(f"Boundary group '{bname}' is empty.")
            gmsh.model.addPhysicalGroup(2, stags, tag=ptag)
            gmsh.model.setPhysicalName(2, ptag, bname)

        gmsh.model.mesh.generate(3)
        # Remove sliver tetrahedra. They set OpenFOAM's Courant time step for the
        # whole mesh; Netgen optimization raises the worst cell quality from about
        # 0.33 to 0.5 and allows roughly 1.5x larger time steps.
        gmsh.model.mesh.optimize("Netgen")
        _verify_every_boundary_has_exactly_one_physical(plasma_vols)
        _verify_all_boundary_faces_mapped(plasma_vols)
        _verify_tetra_only()

        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.write(str(Path(out_msh).resolve()))
    finally:
        gmsh.finalize()

    return boundary_tags


# ─── Structured (graded hexahedral) mesh ───────────────────────────

def graded_coordinates(length, size_start, size_end, size_max, growth_rate):
    """Node coordinates on [0, length] with cell sizes graded from both ends.

    The target size grows geometrically away from each end, by ``growth_rate``
    per cell, from ``size_start`` / ``size_end`` up to ``size_max``. Nodes are
    placed so that the local cell size follows that target; the ends are exact.
    """
    if length <= 0:
        raise ValueError("Segment length must be positive")
    size_start = min(size_start, size_max)
    size_end = min(size_end, size_max)
    slope = growth_rate - 1.0  # geometric growth ~ linear growth of size with distance

    def target(s):
        return min(size_start + slope * s, size_end + slope * (length - s), size_max)

    samples = 4000
    ds = length / samples
    cumulative = [0.0]
    for k in range(samples):
        s_mid = (k + 0.5) * ds
        cumulative.append(cumulative[-1] + ds / target(s_mid))
    # Round up so the configured sizes are upper bounds
    n_cells = max(1, int(math.ceil(cumulative[-1] - 1e-9)))
    scale = n_cells / cumulative[-1]

    coords = [0.0]
    k = 0
    for j in range(1, n_cells):
        goal = j / scale
        while cumulative[k + 1] < goal:
            k += 1
        frac = (goal - cumulative[k]) / (cumulative[k + 1] - cumulative[k])
        coords.append((k + frac) * ds)
    coords.append(length)
    return coords


def _extrusion_heights(coords):
    """Cumulative layer fractions for Gmsh extrusion from node coordinates."""
    length = coords[-1] - coords[0]
    return [(c - coords[0]) / length for c in coords[1:]]


def _verify_hexahedra_only():
    elem_types, _, _ = gmsh.model.mesh.getElements(3)
    bad_types = []
    for et in elem_types:
        name, dim, _, _, _, _ = gmsh.model.mesh.getElementProperties(et)
        if dim != 3 or "hexahedron" not in name.lower():
            bad_types.append((et, name))
    if bad_types:
        raise RuntimeError(
            "3D mesh contains non-hexahedral elements: "
            + ", ".join(f"type={et}:{name}" for et, name in bad_types)
        )


def generate_structured(
    out_msh="channel.msh",
    channel_config=None,
    cell_size=0.0025,
    wall_cell_size=0.00025,
    electrode_edge_cell_size=0.001,
    growth_rate=1.2,
):
    """Generate a graded, structured hexahedral mesh of the channel.

    Cells are finest at the insulating and electrode walls, to resolve the cold
    thermal boundary layer, and streamwise at the electrode edges, where the
    current concentrates; they grow by at most ``growth_rate`` per cell up to
    ``cell_size`` in the core. The channel is extruded along x in segments split
    at every electrode edge, so each electrode is its own boundary surface.
    Boundary names and tags are identical to the tetrahedral generator.
    """
    cfg = {**DEFAULTS, **(channel_config or {})}
    num_pairs = cfg["num_pairs"]
    L = cfg["channel_length"]
    H = cfg["channel_height"]
    W = cfg["channel_width"]
    el = cfg["electrode_length"]
    centers = compute_electrode_centers(L, num_pairs, cfg.get("electrode_centers"))
    tags = build_boundary_tags(num_pairs)

    # Streamwise segments between the inlet, every electrode edge, and the outlet
    tol = 1e-9 * L
    electrode_spans = [(xc - el / 2, xc + el / 2) for xc in centers]
    edges = sorted({0.0, L} | {x for span in electrode_spans for x in span})
    breaks = [edges[0]]
    for x in edges[1:]:
        if x - breaks[-1] > tol:
            breaks.append(x)
    electrode_edges = {x for span in electrode_spans for x in span}

    def is_electrode_edge(x):
        return any(abs(x - e) <= tol for e in electrode_edges)

    def electrode_pair_of_segment(x0, x1):
        for index, (a, b) in enumerate(electrode_spans, start=1):
            if abs(x0 - a) <= tol and abs(x1 - b) <= tol:
                return index
        return None

    y_coords = graded_coordinates(H, wall_cell_size, wall_cell_size, cell_size, growth_rate)
    z_coords = graded_coordinates(W, wall_cell_size, wall_cell_size, cell_size, growth_rate)

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.model.add("linear_hall_channel_structured")
        geo = gmsh.model.geo

        # Cross-section at the inlet (x = 0): point -> y line -> y-z surface
        p0 = geo.addPoint(0.0, 0.0, 0.0)
        line = geo.extrude([(0, p0)], 0, H, 0,
                           numElements=[1] * (len(y_coords) - 1),
                           heights=_extrusion_heights(y_coords), recombine=True)
        line_tag = [t for d, t in line if d == 1][0]
        surf = geo.extrude([(1, line_tag)], 0, 0, W,
                           numElements=[1] * (len(z_coords) - 1),
                           heights=_extrusion_heights(z_coords), recombine=True)
        inlet_surface = [t for d, t in surf if d == 2][0]

        groups = {name: [] for name in tags}
        groups["InletX"].append(inlet_surface)
        volumes = []
        base = inlet_surface
        segments = []
        for x0, x1 in zip(breaks[:-1], breaks[1:]):
            size0 = electrode_edge_cell_size if is_electrode_edge(x0) else cell_size
            size1 = electrode_edge_cell_size if is_electrode_edge(x1) else cell_size
            x_coords = graded_coordinates(x1 - x0, size0, size1, cell_size, growth_rate)
            out = geo.extrude([(2, base)], x1 - x0, 0, 0,
                              numElements=[1] * (len(x_coords) - 1),
                              heights=_extrusion_heights(x_coords), recombine=True)
            top = out[0][1]
            volume = out[1][1]
            laterals = [t for d, t in out[2:] if d == 2]
            volumes.append(volume)
            segments.append((x0, x1, laterals))
            base = top
        groups["OutletX"].append(base)

        geo.synchronize()

        rel = 1e-6 * max(L, H, W)
        for x0, x1, laterals in segments:
            pair = electrode_pair_of_segment(x0, x1)
            for stag in laterals:
                _, ymin, zmin, _, ymax, zmax = gmsh.model.getBoundingBox(2, stag)
                on_y0 = abs(ymin) <= rel and abs(ymax) <= rel
                on_yH = abs(ymin - H) <= rel and abs(ymax - H) <= rel
                if pair is not None and on_y0:
                    groups[f"CathodeSurface_{pair}"].append(stag)
                elif pair is not None and on_yH:
                    groups[f"AnodeSurface_{pair}"].append(stag)
                else:
                    groups["InsulatorSurface"].append(stag)

        gmsh.model.addPhysicalGroup(3, volumes, tag=MATERIAL_TAGS["Plasma"])
        gmsh.model.setPhysicalName(3, MATERIAL_TAGS["Plasma"], "Plasma")
        for name, ptag in tags.items():
            stags = sorted(set(groups[name]))
            if not stags:
                raise RuntimeError(f"Boundary group '{name}' is empty.")
            gmsh.model.addPhysicalGroup(2, stags, tag=ptag)
            gmsh.model.setPhysicalName(2, ptag, name)

        gmsh.model.mesh.generate(3)
        _verify_every_boundary_has_exactly_one_physical(volumes)
        _verify_all_boundary_faces_mapped(volumes)
        _verify_hexahedra_only()

        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.write(str(Path(out_msh).resolve()))
    finally:
        gmsh.finalize()

    return tags


# ─── CLI entry point ───────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate tetrahedral MSH 2.2 mesh for a linear Hall channel."
    )
    parser.add_argument(
        "--out",
        default="channel.msh",
        help="Output mesh path (MSH 2.2 ASCII).",
    )
    parser.add_argument(
        "--mesh-size-min",
        type=float,
        default=None,
        help="Global minimum tetra edge length.",
    )
    parser.add_argument(
        "--mesh-size-max",
        type=float,
        default=None,
        help="Global maximum tetra edge length.",
    )
    parser.add_argument(
        "--mesh-size-factor",
        type=float,
        default=1.0,
        help="Global size scale factor (>1 coarser, <1 finer).",
    )
    args = parser.parse_args()
    generate(
        args.out,
        args.mesh_size_min,
        args.mesh_size_max,
        args.mesh_size_factor,
    )
