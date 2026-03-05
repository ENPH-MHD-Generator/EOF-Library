import argparse
from pathlib import Path

import gmsh


MATERIAL_TAGS = {
    "Plasma": 1,
    "Cathode": 2,
    "Anode": 3,
    "Insulator": 4,
}

BOUNDARY_TAGS = {
    "InletX": 20,
    "OutletX": 21,
    "InsulatorSurface": 30,
    "CathodeSurface": 40,
    "AnodeSurface": 41,
}

STEP_TO_MATERIAL = {
    "air": "Plasma",
    "electrode-": "Cathode",
    "electrode+": "Anode",
    "guide": "Insulator",
}

# Higher rank wins when OCC fragment maps a volume to multiple source materials.
# This is typical when plasma CAD encloses embedded solid inserts.
MATERIAL_PRIORITY = {
    "Plasma": 0,
    "Insulator": 1,
    "Cathode": 2,
    "Anode": 3,
}


def _classify_step_file(step_path: Path) -> str:
    name = step_path.stem.lower()
    for key, material in STEP_TO_MATERIAL.items():
        if key in name:
            return material
    raise RuntimeError(
        f"Cannot map STEP '{step_path.name}' to a material. "
        f"Expected one of: {list(STEP_TO_MATERIAL.keys())}"
    )


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
    # Keep tetrahedral unstructured meshing, but allow global size controls.
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


def _assign_material_physical_volumes(step_dir: Path):
    occ = gmsh.model.occ
    step_files = sorted(step_dir.glob("*.step"))
    if not step_files:
        raise RuntimeError(f"No STEP files found in {step_dir}")

    imported = []
    imported_materials = []
    for step_file in step_files:
        material = _classify_step_file(step_file)
        entities = occ.importShapes(str(step_file))
        vols = [tag for dim, tag in entities if dim == 3]
        if not vols:
            raise RuntimeError(f"STEP file '{step_file.name}' did not create any OCC volumes.")
        for vtag in vols:
            imported.append((3, vtag))
            imported_materials.append(material)

    if not imported:
        raise RuntimeError("No 3D volumes were imported from STEP geometry.")

    # Fragment to enforce conformal interfaces across imported domains.
    # out_map[i] contains the resulting entities generated from imported[i].
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

    # Validate complete and unique coverage of all resulting volumes.
    all_model_vols = {tag for dim, tag in gmsh.model.getEntities(3) if dim == 3}
    vol_to_materials = {v: [] for v in all_model_vols}
    for material, vols in material_vols.items():
        for v in set(vols):
            if v in vol_to_materials:
                vol_to_materials[v].append(material)

    unassigned = sorted(v for v, mats in vol_to_materials.items() if len(mats) == 0)
    if unassigned:
        raise RuntimeError(
            "Post-fragment material mapping is invalid. "
            f"Unassigned volumes: {unassigned}"
        )

    # Resolve multi-mapped volumes deterministically by explicit material priority.
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

    # Rebuild material volumes from resolved one-to-one mapping.
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
    insulator_surfs = []
    anode_surfs = []
    cathode_surfs = []
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
    return (
        sorted(set(insulator_surfs)),
        sorted(set(anode_surfs)),
        sorted(set(cathode_surfs)),
    )


def main(
    step_dir: str,
    out_msh: str,
    mesh_size_min,
    mesh_size_max,
    mesh_size_factor: float,
):
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1)
        _configure_mesh_sizing(mesh_size_min, mesh_size_max, mesh_size_factor)
        gmsh.model.add("channel_step_plasma_only")

        step_path = Path(step_dir).resolve()
        volume_tags, vol_to_material = _assign_material_physical_volumes(step_path)

        # Identify plasma interface surfaces before removing solids.
        ins_surfs, anode_surfs, cathode_surfs = _collect_material_interface_surfaces(vol_to_material)
        if not ins_surfs:
            raise RuntimeError("No Plasma-Insulator interface surfaces found for 'InsulatorSurface'.")
        if not anode_surfs:
            raise RuntimeError("No Plasma-Anode interface surfaces found for 'AnodeSurface'.")
        if not cathode_surfs:
            raise RuntimeError("No Plasma-Cathode interface surfaces found for 'CathodeSurface'.")
        interface_map = {
            "InsulatorSurface": set(ins_surfs),
            "AnodeSurface": set(anode_surfs),
            "CathodeSurface": set(cathode_surfs),
        }

        # Remove solids so the final mesh contains only the plasma region.
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
        for stag in exterior_surfs:
            xmin, _, _, xmax, _, _ = _entity_bbox(2, stag)
            if abs(xmin - xmax) <= tol and abs(xmin - gxmin) <= tol:
                groups["InletX"].append(stag)
                continue
            if abs(xmin - xmax) <= tol and abs(xmax - gxmax) <= tol:
                groups["OutletX"].append(stag)
                continue

            matched = False
            for bname in ("InsulatorSurface", "CathodeSurface", "AnodeSurface"):
                if stag in interface_map[bname]:
                    groups[bname].append(stag)
                    matched = True
                    break
            if not matched:
                raise RuntimeError(
                    "Plasma exterior surface is not inlet/outlet/interface classified: "
                    f"{stag}"
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
        description="Import channel STEP geometry and create tetrahedral MSH 2.2 mesh."
    )
    parser.add_argument(
        "--step-dir",
        default="channel_step",
        help="Directory containing STEP files for the CAD model.",
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
        args.step_dir,
        args.out,
        args.mesh_size_min,
        args.mesh_size_max,
        args.mesh_size_factor,
    )
