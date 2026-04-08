#!/usr/bin/env python3
"""Generate all electrode-dependent files from a centralized electrodes.yaml.

Outputs (relative to --output-dir, default "."):
  - case.sif              Elmer solver input with boundary conditions
  - 0/U, 0/T, 0/p_rgh,   OpenFOAM initial / boundary fields
    0/B, 0/Potential,
    0/J_dens, 0/electric_field, 0/JH
  - channel.msh           Gmsh mesh  (only with --mesh)

Usage:
    python configure.py                       # default electrodes.yaml
    python configure.py --config my.yaml      # custom config
    python configure.py --mesh                # also regenerate the mesh
    python configure.py --dry-run             # preview, don't write
"""

import argparse
import os
import sys
from copy import deepcopy
from pathlib import Path

import yaml


# ===================================================================
# Config loading
# ===================================================================

def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def resolve_pairs(config):
    """Return a list of pair dicts: [{x_center, length, resistance}, ...]."""
    elec = config["electrodes"]
    if "pairs" in elec:
        pairs = elec["pairs"]
        for i, p in enumerate(pairs):
            for key in ("x_center", "length", "resistance"):
                if key not in p:
                    raise ValueError(f"Electrode pair {i+1} missing '{key}'")
        return pairs
    if "count" in elec:
        count = elec["count"]
        length = elec.get("length", 0.010)
        resistance = elec.get("resistance", 1.0)
        ch_len = config["channel"]["length"]
        spacing = ch_len / (count + 1)
        return [
            {"x_center": (i + 1) * spacing, "length": length, "resistance": resistance}
            for i in range(count)
        ]
    raise ValueError(
        "electrodes config must have either 'pairs' (list) or 'count' (int)"
    )


# ===================================================================
# Elmer boundary-index mapping
# ===================================================================

def elmer_boundary_map(num_pairs):
    """Compute the Elmer boundary index for each named boundary.

    ElmerGrid assigns indices in sorted order of the Gmsh physical-group
    tags:  InletX(20)=1, OutletX(21)=2, InsulatorSurface(30)=3,
    then CathodeSurface_i / AnodeSurface_i interleaved from tag 40 upward.
    """
    from generate_mesh import build_boundary_tags
    tags = build_boundary_tags(num_pairs)
    sorted_names = sorted(tags, key=lambda n: tags[n])
    return {name: idx + 1 for idx, name in enumerate(sorted_names)}


# ===================================================================
# case.sif generation
# ===================================================================

_CASE_SIF_HEADER = """\
Header
  CHECK KEYWORDS Warn
  Mesh DB "." "meshElmer"
End

Simulation
  Coordinate System = String "Cartesian 3D"
  Simulation Type = Steady

  Steady State Max Iterations = 1000
  Steady State Min Iterations = 1000

  Output Intervals = 1
  Post File = case.ep
End


Body 1
  Name = "fluid"
  Target Bodies(1) = 1
  Equation = 1
  Material = 1
  Body Force = 1
End

Equation 1
  Name = "EOF_HallChannel"
  Active Solvers(10) = 1 2 3 4 5 6 7 8 9 10
End


! --------------------------------------------------
! 1 Scalar Solvers
! --------------------------------------------------

Solver 1
  Exec Solver = Always
  Equation = "DeclareConductivity"
  Procedure = "AllocateSolver" "AllocateSolver"
  Variable = String "Electric Conductivity"
  Variable DOFs = 1
End

Solver 2
  Exec Solver = Always
  Equation = "DeclareUx"
  Procedure = "AllocateSolver" "AllocateSolver"
  Variable = String "Ux"
  Variable DOFs = 1
End

Solver 3
  Exec Solver = Always
  Equation = "DeclareUy"
  Procedure = "AllocateSolver" "AllocateSolver"
  Variable = String "Uy"
  Variable DOFs = 1
End

Solver 4
  Exec Solver = Always
  Equation = "DeclareUz"
  Procedure = "AllocateSolver" "AllocateSolver"
  Variable = String "Uz"
  Variable DOFs = 1
End

Solver 5
  Exec Solver = Always
  Equation = "DeclareBx"
  Procedure = "AllocateSolver" "AllocateSolver"
  Variable = String "Bx"
  Variable DOFs = 1
End

Solver 6
  Exec Solver = Always
  Equation = "DeclareBy"
  Procedure = "AllocateSolver" "AllocateSolver"
  Variable = String "By"
  Variable DOFs = 1
End

Solver 7
  Exec Solver = Always
  Equation = "DeclareBz"
  Procedure = "AllocateSolver" "AllocateSolver"
  Variable = String "Bz"
  Variable DOFs = 1
End


! --------------------------------------------------
! 2 OpenFOAM -> Elmer
! --------------------------------------------------
Solver 8
  Exec Solver = Always
  Equation = "OpenFOAM2Elmer"
  Procedure = "OpenFOAM2Elmer" "OpenFOAM2ElmerSolver"

  Target Variable 1 = String "Electric Conductivity"
  Target Variable 2 = String "Ux"
  Target Variable 3 = String "Uy"
  Target Variable 4 = String "Uz"
  Target Variable 5 = String "Bx"
  Target Variable 6 = String "By"
  Target Variable 7 = String "Bz"
End


! --------------------------------------------------
! 3 Custom Elmer MHD Solver
! --------------------------------------------------
Solver 9
  Exec Solver = Always
  Equation = "Static Current Solver"
  Procedure = "MHDSolve" "StatCurrentSolver"

  Variable = Potential
  Variable DOFs = 1

  Calculate Volume Current = Logical True
  Calculate Joule Heating  = Logical True

  Nonlinear System Max Iterations = 40
  Nonlinear System Convergence Tolerance = 5.0e-3
  Nonlinear System Relaxation Factor = 0.7
  Nonlinear System Convergence Without Constraints = Logical True

  Linear System Refactorize = Logical True

  Nonlinear System Newton After Iterations = 0
  Nonlinear System Newton After Tolerance  = 1.0e-3

  Linear System Solver = Iterative
  Linear System Iterative Method = GCR
  Linear System GCR Restart = 200
  Linear System Symmetric = False
  Linear System Preconditioning = ILU1
  Linear System Max Iterations = 12000
  Linear System Convergence Tolerance = 1.0e-3
  Linear System Abort Not Converged = False
  Linear System Scaling = False
  Linear System Residual Output = 500
End



! --------------------------------------------------
! 4 Elmer -> OpenFOAM: export computed results
! --------------------------------------------------
Solver 10
  Exec Solver = Always
  Equation = "Elmer2OpenFOAM"
  Procedure = "Elmer2OpenFOAM" "Elmer2OpenFOAMSolver"

  Target Variable 1 = String "Volume Current 1"
  Target Variable 2 = String "Volume Current 2"
  Target Variable 3 = String "Volume Current 3"
  Target Variable 4 = String "Joule Heating"
  Target Variable 5 = String "Potential"
End


Material 1
  Name = "Conducting Fluid"

  Electric Conductivity = Variable "Electric Conductivity"
    Real MATC "tx"
End


Body Force 1
  Name = "NoSource"
End
"""


def _generate_case_sif(pairs, bnd_map):
    """Return the full case.sif content as a string."""
    num_pairs = len(pairs)
    lines = [_CASE_SIF_HEADER]

    # Boundary-name -> Elmer index comment block
    lines.append("! -------------------------")
    lines.append("! Boundary Conditions")
    lines.append("! -------------------------")
    lines.append("!   Boundary name            Elmer index")
    for name in sorted(bnd_map, key=lambda n: bnd_map[n]):
        lines.append(f"!   {name:<26s} {bnd_map[name]}")
    lines.append("")

    bc_num = 0

    # Electrode pairs
    for i, pair in enumerate(pairs, start=1):
        bc_num += 1
        lines.append(f"Boundary Condition {bc_num}")
        lines.append(f"  ! CathodeSurface_{i}")
        lines.append(f"  Target Boundaries(1) = {bnd_map[f'CathodeSurface_{i}']}")
        lines.append(f"  Electrode Pair = Integer {i}")
        lines.append(f'  Electrode Sign = String "minus"')
        lines.append(f"  Electrode Resistance = Real {pair['resistance']}")
        lines.append("End")
        lines.append("")

        bc_num += 1
        lines.append(f"Boundary Condition {bc_num}")
        lines.append(f"  ! AnodeSurface_{i}")
        lines.append(f"  Target Boundaries(1) = {bnd_map[f'AnodeSurface_{i}']}")
        lines.append(f"  Electrode Pair = Integer {i}")
        lines.append(f'  Electrode Sign = String "plus"')
        lines.append(f"  Electrode Resistance = Real {pair['resistance']}")
        lines.append("End")
        lines.append("")

    # Insulator
    bc_num += 1
    lines.append(f"Boundary Condition {bc_num}")
    lines.append("  ! InsulatorSurface")
    lines.append(f"  Target Boundaries(1) = {bnd_map['InsulatorSurface']}")
    lines.append("End")
    lines.append("")

    # Inlet
    bc_num += 1
    lines.append(f"Boundary Condition {bc_num}")
    lines.append("  ! InletX")
    lines.append(f"  Target Boundaries(1) = {bnd_map['InletX']}")
    lines.append("End")
    lines.append("")

    # Outlet
    bc_num += 1
    lines.append(f"Boundary Condition {bc_num}")
    lines.append("  ! OutletX")
    lines.append(f"  Target Boundaries(1) = {bnd_map['OutletX']}")
    lines.append("End")
    lines.append("")

    return "\n".join(lines)


# ===================================================================
# OpenFOAM 0/ file generation
# ===================================================================

_OF_HEADER = """\
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
    class       {field_class};
    location    "0";
    object      {field_name};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

{comment}dimensions {dimensions};

internalField   {internal_field};


boundaryField
{{
{patches}}}


// ************************************************************************* //
"""


def _fmt_vec(v):
    """Format a list/tuple as an OpenFOAM vector string: (x y z)."""
    return "(" + " ".join(str(x) for x in v) + ")"


def _render_patch(name, bc_dict):
    """Render one boundaryField patch entry."""
    lines = [f"    {name}", "    {"]
    lines.append(f"        type    {bc_dict['type']};")
    if "value" in bc_dict:
        lines.append(f"        value   {bc_dict['value']};")
    lines.append("    }")
    return "\n".join(lines)


def _resolve_bc(bc_dict, fmt):
    """Return a copy of bc_dict with format placeholders resolved."""
    resolved = {}
    for k, v in bc_dict.items():
        resolved[k] = v.format_map(fmt) if isinstance(v, str) else v
    return resolved


def _generate_openfoam_field(field_def, electrode_patches, fmt):
    """Generate one OpenFOAM 0/ file content."""
    patches_block = []

    patches_block.append(_render_patch("InletX", _resolve_bc(field_def["inlet"], fmt)))
    patches_block.append("")
    patches_block.append(_render_patch("OutletX", _resolve_bc(field_def["outlet"], fmt)))
    patches_block.append("")
    patches_block.append(_render_patch("InsulatorSurface", _resolve_bc(field_def["wall"], fmt)))

    for pname in electrode_patches:
        patches_block.append("")
        patches_block.append(_render_patch(pname, _resolve_bc(field_def["wall"], fmt)))

    patches_block.append("")
    patches_block.append(_render_patch("defaultFaces", _resolve_bc(field_def["wall"], fmt)))
    patches_block.append("")

    comment_line = ""
    if field_def.get("comment"):
        comment_line = field_def["comment"] + "\n"

    internal = field_def["internal_field"].format_map(fmt)

    return _OF_HEADER.format(
        field_class=field_def["class"],
        field_name=field_def["name"],
        comment=comment_line,
        dimensions=field_def["dimensions"],
        internal_field=internal,
        patches="\n".join(patches_block),
    )


def _openfoam_field_defs():
    """Return the list of OpenFOAM field definitions.

    Placeholders like {inlet_velocity} are resolved at generation time via
    a format-mapping dict built from the physics config.
    """
    return [
        {
            "name": "U",
            "class": "volVectorField",
            "dimensions": "[0 1 -1 0 0 0 0]",
            "internal_field": "uniform (0 0 0)",
            "comment": None,
            "inlet": {"type": "fixedValue", "value": "uniform {inlet_velocity}"},
            "outlet": {"type": "zeroGradient"},
            "wall": {"type": "noSlip"},
        },
        {
            "name": "T",
            "class": "volScalarField",
            "dimensions": "[0 0 0 1 0 0 0]",
            "internal_field": "uniform 100",
            "comment": "// Temperature [K]",
            "inlet": {"type": "fixedValue", "value": "uniform {inlet_temperature}"},
            "outlet": {"type": "zeroGradient"},
            "wall": {"type": "zeroGradient"},
        },
        {
            "name": "p_rgh",
            "class": "volScalarField",
            "dimensions": "[1 -1 -2 0 0 0 0]",
            "internal_field": "uniform 0",
            "comment": "// Dynamic pressure (p - rho*g*h)",
            "inlet": {"type": "zeroGradient"},
            "outlet": {"type": "fixedValue", "value": "uniform 0"},
            "wall": {"type": "zeroGradient"},
        },
        {
            "name": "B",
            "class": "volVectorField",
            "dimensions": "[1 0 -2 0 0 -1 0]",
            "internal_field": "uniform {B_field}",
            "comment": "// Magnetic flux density [Tesla]",
            "inlet": {"type": "zeroGradient"},
            "outlet": {"type": "zeroGradient"},
            "wall": {"type": "zeroGradient"},
        },
        {
            "name": "Potential",
            "class": "volScalarField",
            "dimensions": "[1 2 -3 0 0 -1 0]",
            "internal_field": "uniform 0",
            "comment": "// Electric potential [V]",
            "inlet": {"type": "zeroGradient"},
            "outlet": {"type": "zeroGradient"},
            "wall": {"type": "zeroGradient"},
        },
        {
            "name": "J_dens",
            "class": "volVectorField",
            "dimensions": "[0 -2 0 0 0 1 0]",
            "internal_field": "uniform (0 0 0)",
            "comment": "// Volume current density [A/m^2]",
            "inlet": {"type": "zeroGradient"},
            "outlet": {"type": "zeroGradient"},
            "wall": {"type": "zeroGradient"},
        },
        {
            "name": "electric_field",
            "class": "volVectorField",
            "dimensions": "[1 1 -3 0 0 -1 0]",
            "internal_field": "uniform (0 0 0)",
            "comment": "// Electric field E = -grad(Potential) [V/m]",
            "inlet": {"type": "zeroGradient"},
            "outlet": {"type": "zeroGradient"},
            "wall": {"type": "zeroGradient"},
        },
        {
            "name": "JH",
            "class": "volScalarField",
            "dimensions": "[1 -1 -3 0 0 0 0]",
            "internal_field": "uniform 0",
            "comment": "// Joule heating power density [W/m^3]",
            "inlet": {"type": "zeroGradient"},
            "outlet": {"type": "zeroGradient"},
            "wall": {"type": "zeroGradient"},
        },
    ]


# ===================================================================
# Mesh generation (delegates to generate_mesh.py)
# ===================================================================

def _generate_mesh(config, pairs, output_dir):
    from generate_mesh import generate as gen_mesh

    # All pairs currently share the electrode_length of the first pair.
    # Per-pair lengths would require generate_mesh changes.
    electrode_length = pairs[0]["length"]
    centers = [p["x_center"] for p in pairs]

    ch = config["channel"]
    mesh_cfg = config.get("mesh", {})

    channel_config = {
        "num_pairs": len(pairs),
        "channel_length": ch["length"],
        "channel_height": ch["height"],
        "channel_width": ch["width"],
        "electrode_length": electrode_length,
        "wall_thickness": ch["wall_thickness"],
        "electrode_centers": centers,
    }

    out_msh = str(Path(output_dir) / "channel.msh")
    gen_mesh(
        out_msh=out_msh,
        mesh_size_min=mesh_cfg.get("size_min"),
        mesh_size_max=mesh_cfg.get("size_max"),
        mesh_size_factor=mesh_cfg.get("size_factor", 1.0),
        channel_config=channel_config,
    )
    print(f"  wrote {out_msh}")


# ===================================================================
# Orchestrator
# ===================================================================

def configure(config, output_dir=".", do_mesh=False, dry_run=False):
    """Generate all files from the given config dict."""
    pairs = resolve_pairs(config)
    num_pairs = len(pairs)
    physics = config.get("physics", {})

    print(f"Electrode pairs: {num_pairs}")
    for i, p in enumerate(pairs, 1):
        print(f"  {i}: x={p['x_center']:.4f} m, L={p['length']:.4f} m, R={p['resistance']} Ω")

    # ── Format-mapping for OpenFOAM placeholders ──
    fmt = {
        "inlet_velocity": _fmt_vec(physics.get("inlet_velocity", [0, 0, 0])),
        "inlet_temperature": physics.get("inlet_temperature", 300),
        "B_field": _fmt_vec(physics.get("B_field", [0, 0, 0])),
    }

    # ── Ordered electrode patch names ──
    electrode_patches = []
    for i in range(1, num_pairs + 1):
        electrode_patches.append(f"CathodeSurface_{i}")
        electrode_patches.append(f"AnodeSurface_{i}")

    # ── Elmer boundary index map ──
    bnd_map = elmer_boundary_map(num_pairs)

    # ── case.sif ──
    sif_content = _generate_case_sif(pairs, bnd_map)
    sif_path = Path(output_dir) / "case.sif"
    if dry_run:
        print(f"\n[dry-run] would write {sif_path} ({len(sif_content)} chars)")
    else:
        sif_path.write_text(sif_content)
        print(f"  wrote {sif_path}")

    # ── OpenFOAM 0/ files ──
    zero_dir = Path(output_dir) / "0"
    zero_dir.mkdir(exist_ok=True)
    for field_def in _openfoam_field_defs():
        content = _generate_openfoam_field(field_def, electrode_patches, fmt)
        fpath = zero_dir / field_def["name"]
        if dry_run:
            print(f"[dry-run] would write {fpath}")
        else:
            fpath.write_text(content)
            print(f"  wrote {fpath}")

    # ── Mesh ──
    if do_mesh:
        if dry_run:
            print("[dry-run] would regenerate mesh")
        else:
            _generate_mesh(config, pairs, output_dir)

    print("Done.")


# ===================================================================
# CLI
# ===================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate case.sif, OpenFOAM 0/ files, and optionally the mesh "
        "from a centralized electrode configuration.",
    )
    parser.add_argument(
        "--config",
        default="electrodes.yaml",
        help="Path to the YAML config file (default: electrodes.yaml).",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory to write generated files into (default: cwd).",
    )
    parser.add_argument(
        "--mesh",
        action="store_true",
        help="Also regenerate the Gmsh mesh (channel.msh).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be generated without writing files.",
    )
    args = parser.parse_args()

    config = load_config(args.config)
    configure(config, args.output_dir, args.mesh, args.dry_run)


if __name__ == "__main__":
    main()
