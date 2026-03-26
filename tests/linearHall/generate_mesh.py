import argparse
from pathlib import Path

import gmsh


# ─── Geometry constants (edit these to change the channel layout) ─
NUM_ELECTRODE_PAIRS = 4
CHANNEL_LENGTH      = 0.200    # m  (200 mm)
CHANNEL_HEIGHT      = 0.050    # m  (10 mm)
CHANNEL_WIDTH       = 0.050    # m  (10 mm)
ELECTRODE_LENGTH    = 0.010    # m  (10 mm along flow direction)
WALL_THICKNESS      = 0.005    # m  (2 mm, used for insulator / electrode shells)

# ─── Material tags (3-D physical groups) ──────────────────────────
MATERIAL_TAGS = {
    "Plasma": 1,
    "Cathode": 2,
    "Anode": 3,
    "Insulator": 4,
}

# Higher rank wins when OCC fragment maps a volume to multiple source materials.
MATERIAL_PRIORITY = {
    "Plasma": 0,
    "Insulator": 1,
    "Cathode": 2,
    "Anode": 3,
}

# ─── Boundary tags (2-D physical groups) ──────────────────────────

def _build_boundary_tags():
    tags = {
        "InletX": 20,
        "OutletX": 21,
        "InsulatorSurface": 30,
    }
    for i in range(NUM_ELECTRODE_PAIRS):
        tags[f"CathodeSurface_{i + 1}"] = 40 + 2 * i
        tags[f"AnodeSurface_{i + 1}"] = 41 + 2 * i
    return tags

BOUNDARY_TAGS = _build_boundary_tags()


def _electrode_centers():
    spacing = CHANNEL_LENGTH / (NUM_ELECTRODE_PAIRS + 1)
    return [(i + 1) * spacing for i in range(NUM_ELECTRODE_PAIRS)]


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
    # Exterior surface belongs to exactly one volume.
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


def _build_geometry_and_assign_materials():
    """Build a rectangular Hall channel with evenly-spaced electrode pairs."""
    occ = gmsh.model.occ
    L, H, W, t = CHANNEL_LENGTH, CHANNEL_HEIGHT, CHANNEL_WIDTH, WALL_THICKNESS
    el = ELECTRODE_LENGTH

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

    centers = _electrode_centers()
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


def _collect_material_interface_surfaces(vol_to_material):
    """Find interface surfaces between plasma and surrounding materials.

    Returns (insulator_surfs, cathode_groups, anode_groups) where the
    electrode groups are dicts mapping pair index (0-based) to surface tags.
    """
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

    centers = _electrode_centers()
    cathode_groups = {i: [] for i in range(NUM_ELECTRODE_PAIRS)}
    anode_groups = {i: [] for i in range(NUM_ELECTRODE_PAIRS)}

    for stag in cathode_surfs:
        xmin, _, _, xmax, _, _ = gmsh.model.getBoundingBox(2, stag)
        x_mid = (xmin + xmax) / 2
        pair_idx = min(range(NUM_ELECTRODE_PAIRS),
                       key=lambda j: abs(x_mid - centers[j]))
        cathode_groups[pair_idx].append(stag)

    for stag in anode_surfs:
        xmin, _, _, xmax, _, _ = gmsh.model.getBoundingBox(2, stag)
        x_mid = (xmin + xmax) / 2
        pair_idx = min(range(NUM_ELECTRODE_PAIRS),
                       key=lambda j: abs(x_mid - centers[j]))
        anode_groups[pair_idx].append(stag)

    return sorted(set(insulator_surfs)), cathode_groups, anode_groups


def main(
    out_msh: str,
    mesh_size_min,
    mesh_size_max,
    mesh_size_factor: float,
):
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1)
        _configure_mesh_sizing(mesh_size_min, mesh_size_max, mesh_size_factor)
        gmsh.model.add("linear_hall_channel")

        volume_tags, vol_to_material = _build_geometry_and_assign_materials()

        ins_surfs, cathode_groups, anode_groups = _collect_material_interface_surfaces(
            vol_to_material
        )
        if not ins_surfs:
            raise RuntimeError("No Plasma-Insulator interface surfaces found.")

        interface_set = set(ins_surfs)
        for i in range(NUM_ELECTRODE_PAIRS):
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

        groups = {name: [] for name in BOUNDARY_TAGS}
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
            for i in range(NUM_ELECTRODE_PAIRS):
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

        for bname, ptag in BOUNDARY_TAGS.items():
            stags = sorted(set(groups[bname]))
            if not stags:
                raise RuntimeError(f"Boundary group '{bname}' is empty.")
            gmsh.model.addPhysicalGroup(2, stags, tag=ptag)
            gmsh.model.setPhysicalName(2, ptag, bname)

        gmsh.model.mesh.generate(3)
        _verify_every_boundary_has_exactly_one_physical(plasma_vols)
        _verify_all_boundary_faces_mapped(plasma_vols)
        _verify_tetra_only()

        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.write(str(Path(out_msh).resolve()))
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate tetrahedral MSH 2.2 mesh for a linear Hall channel "
        f"with {NUM_ELECTRODE_PAIRS} electrode pairs."
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
        help="Global minimum tetra edge length (smaller -> finer mesh).",
    )
    parser.add_argument(
        "--mesh-size-max",
        type=float,
        default=None,
        help="Global maximum tetra edge length (smaller -> finer mesh).",
    )
    parser.add_argument(
        "--mesh-size-factor",
        type=float,
        default=1.0,
        help="Global size scale factor (>1 coarser, <1 finer).",
    )
    args = parser.parse_args()
    main(
        args.out,
        args.mesh_size_min,
        args.mesh_size_max,
        args.mesh_size_factor,
    )
