import argparse
from collections import defaultdict
from pathlib import Path
from typing import Optional, Tuple

import gmsh


MATERIAL_TAGS = {
    "Plasma": 1,
    "Cathode": 2,
    "Anode": 3,
    "Insulator": 4,
}

# Physical tags (2D). ElmerGrid maps boundaries by increasing tag -> BC index 1..N.
BOUNDARY_TAGS = {
    "InletX": 20,
    "OutletX": 21,
    "InsulatorSurface": 30,
    "CathodeSurface_1": 40,
    "CathodeSurface_2": 41,
    "CathodeSurface_3": 42,
    "AnodeSurface_1": 43,
    "AnodeSurface_2": 44,
    "AnodeSurface_3": 45,
}

ELECTRODE_SLOTS = (1, 2, 3)

STEP_TO_MATERIAL = {
    "air": "Plasma",
    "plasma": "Plasma",
    "fluid": "Plasma",
    "gas": "Plasma",
    "cathode_1": "Cathode",
    "cathode_2": "Cathode",
    "cathode_3": "Cathode",
    "anode_1": "Anode",
    "anode_2": "Anode",
    "anode_3": "Anode",
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


def _classify_step_file(step_path: Path) -> Tuple[str, Optional[int]]:
    """Map STEP stem to (material, electrode slot 1..3 or None)."""
    name = step_path.stem.lower()
    keys = sorted(STEP_TO_MATERIAL.keys(), key=lambda k: (-len(k), k))
    for key in keys:
        i = name.find(key)
        if i < 0:
            continue
        j = i + len(key)
        if j < len(name) and name[j].isdigit():
            continue
        material = STEP_TO_MATERIAL[key]
        slot: Optional[int] = None
        if material in ("Cathode", "Anode"):
            tail = key.rsplit("_", 1)[-1]
            if not tail.isdigit():
                raise RuntimeError(
                    f"Electrode STEP key '{key}' must end with _<n> (e.g. cathode_2)."
                )
            slot = int(tail)
        return material, slot
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


def _assign_material_physical_volumes(step_dir: Path, *, auto_plasma_box: bool):
    occ = gmsh.model.occ
    step_files = sorted(step_dir.glob("*.step"))
    if not step_files:
        raise RuntimeError(f"No STEP files found in {step_dir}")

    imported = []
    imported_meta: list[Tuple[str, Optional[int]]] = []
    for step_file in step_files:
        material, slot = _classify_step_file(step_file)
        entities = occ.importShapes(str(step_file))
        vols = [tag for dim, tag in entities if dim == 3]
        if not vols:
            raise RuntimeError(f"STEP file '{step_file.name}' did not create any OCC volumes.")
        for vtag in vols:
            imported.append((3, vtag))
            imported_meta.append((material, slot))

    if not imported:
        raise RuntimeError("No 3D volumes were imported from STEP geometry.")

    occ.synchronize()
    if not any(m == "Plasma" for m, _ in imported_meta):
        if not auto_plasma_box:
            raise RuntimeError(
                "No STEP maps to Plasma (filename must contain one of: "
                f"{[k for k, v in STEP_TO_MATERIAL.items() if v == 'Plasma']}). "
                "Export the gas volume from CAD, or run with default "
                "auto-plasma-box enabled (omit --no-auto-plasma-box)."
            )
        print(
            "[gmsh] WARNING: No air/plasma STEP; using a tight bounding box as the fluid "
            "domain. Prefer exporting the cavity as a STEP with 'plasma' or 'air' in the name."
        )
        gxmin, gymin, gzmin, gxmax, gymax, gzmax = gmsh.model.getBoundingBox(-1, -1)
        dom = max(gxmax - gxmin, gymax - gymin, gzmax - gzmin)
        pad = max(1e-9, 1e-6 * dom)
        x0, y0, z0 = gxmin - pad, gymin - pad, gzmin - pad
        dx = (gxmax - gxmin) + 2 * pad
        dy = (gymax - gymin) + 2 * pad
        dz = (gzmax - gzmin) + 2 * pad
        box_tag = occ.addBox(x0, y0, z0, dx, dy, dz)
        occ.synchronize()
        imported.insert(0, (3, box_tag))
        imported_meta.insert(0, ("Plasma", None))

    # Fragment to enforce conformal interfaces across imported domains.
    # out_map[i] contains the resulting entities generated from imported[i].
    _, out_map = occ.fragment(imported, [])
    occ.synchronize()

    material_vols = {name: [] for name in MATERIAL_TAGS}
    if len(out_map) != len(imported):
        raise RuntimeError(
            f"Unexpected OCC fragment map size: got {len(out_map)}, expected {len(imported)}."
        )

    vol_claims: dict[int, list[Tuple[str, Optional[int]]]] = defaultdict(list)
    for i, children in enumerate(out_map):
        material, slot = imported_meta[i]
        for dim, tag in children:
            if dim == 3:
                material_vols[material].append(tag)
                vol_claims[tag].append((material, slot))

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

    vol_to_electrode_slot: dict[int, Optional[int]] = {}
    for v in all_model_vols:
        mat = resolved_vol_to_material[v]
        claims = vol_claims.get(v, [])
        if mat in ("Cathode", "Anode"):
            slots = sorted({sl for m, sl in claims if m == mat and sl is not None})
            if len(slots) != 1:
                raise RuntimeError(
                    f"Volume {v} resolved as {mat} but expected exactly one electrode index "
                    f"from fragment claims {claims}; got slots {slots}"
                )
            vol_to_electrode_slot[v] = slots[0]
        else:
            vol_to_electrode_slot[v] = None

    new_vols = sorted(all_model_vols)
    if not new_vols:
        raise RuntimeError("No volumes present after OCC fragment/synchronize.")
    return sorted(new_vols), vol_to_material, vol_to_electrode_slot


def _collect_electrode_interface_surfaces(vol_to_material, vol_to_electrode_slot):
    insulator_surfs = []
    cathode_by_slot = {s: [] for s in ELECTRODE_SLOTS}
    anode_by_slot = {s: [] for s in ELECTRODE_SLOTS}
    for _, stag in gmsh.model.getEntities(2):
        up, _ = gmsh.model.getAdjacencies(2, stag)
        if len(up) != 2:
            continue
        m0, m1 = vol_to_material.get(up[0]), vol_to_material.get(up[1])
        mats = {m0, m1}
        if mats == {"Plasma", "Insulator"}:
            insulator_surfs.append(stag)
        elif mats == {"Plasma", "Cathode"}:
            cv = up[0] if m0 == "Cathode" else up[1]
            slot = vol_to_electrode_slot[cv]
            if slot is None:
                raise RuntimeError(f"Plasma-Cathode face {stag} has no electrode slot for volume {cv}")
            cathode_by_slot[slot].append(stag)
        elif mats == {"Plasma", "Anode"}:
            av = up[0] if m0 == "Anode" else up[1]
            slot = vol_to_electrode_slot[av]
            if slot is None:
                raise RuntimeError(f"Plasma-Anode face {stag} has no electrode slot for volume {av}")
            anode_by_slot[slot].append(stag)
    cathode_by_slot = {s: sorted(set(v)) for s, v in cathode_by_slot.items()}
    anode_by_slot = {s: sorted(set(v)) for s, v in anode_by_slot.items()}
    return sorted(set(insulator_surfs)), cathode_by_slot, anode_by_slot


def main(
    step_dir: str,
    out_msh: str,
    mesh_size_min,
    mesh_size_max,
    mesh_size_factor: float,
    *,
    auto_plasma_box: bool = True,
):
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1)
        _configure_mesh_sizing(mesh_size_min, mesh_size_max, mesh_size_factor)
        gmsh.model.add("channel_step_plasma_only")

        step_path = Path(step_dir).resolve()
        volume_tags, vol_to_material, vol_to_electrode_slot = _assign_material_physical_volumes(
            step_path, auto_plasma_box=auto_plasma_box
        )

        # Identify plasma interface surfaces before removing solids (per electrode index).
        ins_surfs, cathode_by_slot, anode_by_slot = _collect_electrode_interface_surfaces(
            vol_to_material, vol_to_electrode_slot
        )
        if not ins_surfs:
            raise RuntimeError("No Plasma-Insulator interface surfaces found for 'InsulatorSurface'.")
        for s in ELECTRODE_SLOTS:
            if not anode_by_slot[s]:
                raise RuntimeError(
                    f"No Plasma-Anode interface surfaces found for 'AnodeSurface_{s}'."
                )
            if not cathode_by_slot[s]:
                raise RuntimeError(
                    f"No Plasma-Cathode interface surfaces found for 'CathodeSurface_{s}'."
                )
        interface_map = {"InsulatorSurface": set(ins_surfs)}
        for s in ELECTRODE_SLOTS:
            interface_map[f"CathodeSurface_{s}"] = set(cathode_by_slot[s])
            interface_map[f"AnodeSurface_{s}"] = set(anode_by_slot[s])

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
        iface_check_order = ["InsulatorSurface"] + [
            f"CathodeSurface_{s}" for s in ELECTRODE_SLOTS
        ] + [f"AnodeSurface_{s}" for s in ELECTRODE_SLOTS]
        fallback_insulator = []
        for stag in exterior_surfs:
            xmin, _, _, xmax, _, _ = _entity_bbox(2, stag)
            if abs(xmin - xmax) <= tol and abs(xmin - gxmin) <= tol:
                groups["InletX"].append(stag)
                continue
            if abs(xmin - xmax) <= tol and abs(xmax - gxmax) <= tol:
                groups["OutletX"].append(stag)
                continue

            matched = False
            for bname in iface_check_order:
                if stag in interface_map[bname]:
                    groups[bname].append(stag)
                    matched = True
                    break
            if not matched:
                # e.g. y/z faces of the auto-plasma bounding box — treat as wall (insulator BC).
                groups["InsulatorSurface"].append(stag)
                fallback_insulator.append(stag)

        if fallback_insulator:
            print(
                "[gmsh] WARNING: "
                f"{len(fallback_insulator)} plasma boundary face(s) tagged InsulatorSurface "
                "(not inlet/outlet/electrode); typical for auto-plasma-box side caps."
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
        default="channel_step_taper",
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
    parser.add_argument(
        "--no-auto-plasma-box",
        action="store_true",
        help="Fail if no STEP maps to Plasma instead of filling with a bounding box.",
    )
    args = parser.parse_args()
    main(
        args.step_dir,
        args.out,
        args.mesh_size_min,
        args.mesh_size_max,
        args.mesh_size_factor,
        auto_plasma_box=not args.no_auto_plasma_box,
    )
