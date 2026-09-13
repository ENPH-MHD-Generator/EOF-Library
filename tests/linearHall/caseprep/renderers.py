"""Render validated case models into Elmer and OpenFOAM input files."""

from __future__ import annotations

from dataclasses import asdict
from typing import Dict, Iterable, List, Mapping

from .boundaries import INLET, INSULATOR, OUTLET, electrode_patch_names
from .config import CaseConfig


CASE_SIF_HEADER = """\
Header
  CHECK KEYWORDS Warn
  Mesh DB "." "meshElmer"
End

Simulation
  Coordinate System = String "Cartesian 3D"
  Simulation Type = Steady

  ! One iteration per OpenFOAM coupling update; OpenFOAM ends the run by
  ! sending a final status, so these must exceed the number of updates
  Steady State Max Iterations = 100000000
  Steady State Min Iterations = 100000000

  ! No Elmer result files: results are viewed through OpenFOAM, and an ElmerPost
  ! file here would be rewritten every coupling update (tens of MB per run)
  Output Intervals = 0
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
  Active Solvers(12) = 1 2 3 4 5 6 7 8 9 10 11 12
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

! Pressure and gas temperature feed the seeded-plasma (Saha) model

Solver 8
  Exec Solver = Always
  Equation = "DeclarePressure"
  Procedure = "AllocateSolver" "AllocateSolver"
  Variable = String "Pressure"
  Variable DOFs = 1
End

Solver 9
  Exec Solver = Always
  Equation = "DeclareGasTemp"
  Procedure = "AllocateSolver" "AllocateSolver"
  Variable = String "Gas Temperature"
  Variable DOFs = 1
End

! --------------------------------------------------
! 2 OpenFOAM -> Elmer
! --------------------------------------------------
Solver 10
  Exec Solver = Always
  Equation = "OpenFOAM2Elmer"
  Procedure = "OpenFOAM2Elmer" "OpenFOAM2ElmerSolver"

  Target Variable 1 = String "Ux"
  Target Variable 2 = String "Uy"
  Target Variable 3 = String "Uz"
  Target Variable 4 = String "Bx"
  Target Variable 5 = String "By"
  Target Variable 6 = String "Bz"
  Target Variable 7 = String "Pressure"
  Target Variable 8 = String "Gas Temperature"
End


! --------------------------------------------------
! 3 Custom Elmer MHD Solver
! --------------------------------------------------
Solver 11
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
Solver 12
  Exec Solver = Always
  Equation = "Elmer2OpenFOAM"
  Procedure = "Elmer2OpenFOAM" "Elmer2OpenFOAMSolver"

  ! Element-wise: nodal averaging would smear the electrode-edge peaks of the
  ! current and heating, inflating what OpenFOAM integrates
  Target Variable 1 = String "Element Volume Current 1"
  Target Variable 2 = String "Element Volume Current 2"
  Target Variable 3 = String "Element Volume Current 3"
  Target Variable 4 = String "Element Joule Heating"
  Target Variable 5 = String "Potential"
  Target Variable 6 = String "Electric Conductivity"
  Target Variable 7 = String "Ionization Fraction"
  Target Variable 8 = String "Electron Temperature"
End


Body Force 1
  Name = "NoSource"
End
"""


MATERIAL_TEMPLATE = """\
Material 1
  Name = "Conducting Fluid"

  ! Computed by MHDSolve from the Saha equation for the seed species
  Electric Conductivity = Variable "Electric Conductivity"
    Real MATC "tx"

  Seed Mole Fraction = Real {seed_mole_fraction:g}
  Seed Ionization Energy = Real {seed_ionization_energy:g}  ! eV
  Seed Statistical Weight Ratio = Real {seed_gi_over_gn:g}  ! g_ion / g_neutral
  Seed Electron Neutral Cross Section = Real {seed_cross_section:g}  ! m^2
  Carrier Electron Neutral Cross Section = Real {carrier_cross_section:g}  ! m^2
  Reference Pressure = Real {reference_pressure:g}  ! Pa, added to OpenFOAM gauge p

  Sigma Min = Real {sigma_min:g}
  Sigma Max = Real {sigma_max:g}

  ! Electron temperature: Joule heating vs. elastic losses (Kerrebrock),
  ! otherwise Te = Tgas
  Two Temperature = Logical {two_temperature}
  Carrier Molar Mass = Real {carrier_molar_mass:g}  ! g/mol
  Seed Molar Mass = Real {seed_molar_mass:g}  ! g/mol
  Electron Energy Loss Factor = Real {energy_loss_factor:g}
  Electron Temperature Max = Real {electron_temperature_max:g}  ! K
  Electron Temperature Relaxation = Real {electron_temperature_relaxation:g}
End
"""


class ElmerCaseRenderer:
    """Render an Elmer SIF from a resolved case model."""

    def render(self, config: CaseConfig, boundary_indices: Mapping[str, int]) -> str:
        material = asdict(config.plasma)
        material["two_temperature"] = "True" if config.plasma.two_temperature else "False"
        lines = [CASE_SIF_HEADER, MATERIAL_TEMPLATE.format_map(material)]
        lines.extend(
            (
                "! -------------------------",
                "! Boundary Conditions",
                "! -------------------------",
                "!   Boundary name            Elmer index",
            )
        )
        for name in sorted(boundary_indices, key=boundary_indices.get):
            lines.append(f"!   {name:<26s} {boundary_indices[name]}")
        lines.append("")

        condition_number = 0
        for pair_number, pair in enumerate(config.electrodes.pairs, start=1):
            for role, sign in (("Cathode", "minus"), ("Anode", "plus")):
                condition_number += 1
                name = f"{role}Surface_{pair_number}"
                lines.extend(
                    (
                        f"Boundary Condition {condition_number}",
                        f"  ! {name}",
                        f"  Target Boundaries(1) = {boundary_indices[name]}",
                        f"  Electrode Pair = Integer {pair_number}",
                        f'  Electrode Sign = String "{sign}"',
                        f"  Electrode Resistance = Real {pair.resistance}",
                        "End",
                        "",
                    )
                )

        for name in (INSULATOR, INLET, OUTLET):
            condition_number += 1
            lines.extend(
                (
                    f"Boundary Condition {condition_number}",
                    f"  ! {name}",
                    f"  Target Boundaries(1) = {boundary_indices[name]}",
                    "End",
                    "",
                )
            )
        return "\n".join(lines)


OPENFOAM_HEADER = """\
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


FIELD_DEFINITIONS = (
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
        "internal_field": "uniform {inlet_temperature}",
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
    # Plasma state computed by Elmer. The solver does not read these; they
    # exist so the fields are present in the t = 0 output for viewing.
    {
        "name": "Te",
        "class": "volScalarField",
        "dimensions": "[0 0 0 1 0 0 0]",
        "internal_field": "uniform {inlet_temperature}",
        "comment": "// Electron temperature [K]",
        "inlet": {"type": "zeroGradient"},
        "outlet": {"type": "zeroGradient"},
        "wall": {"type": "zeroGradient"},
    },
    {
        "name": "ionizationFraction",
        "class": "volScalarField",
        "dimensions": "[0 0 0 0 0 0 0]",
        "internal_field": "uniform 0",
        "comment": "// Electron mole fraction n_e / n_heavy",
        "inlet": {"type": "zeroGradient"},
        "outlet": {"type": "zeroGradient"},
        "wall": {"type": "zeroGradient"},
    },
    {
        "name": "elcond_elmer",
        "class": "volScalarField",
        "dimensions": "[-1 -3 3 0 0 2 0]",
        "internal_field": "uniform 0",
        "comment": "// Electrical conductivity from the Saha model [S/m]",
        "inlet": {"type": "zeroGradient"},
        "outlet": {"type": "zeroGradient"},
        "wall": {"type": "zeroGradient"},
    },
)


def _format_vector(values: Iterable[float]) -> str:
    return "(" + " ".join(f"{value:g}" for value in values) + ")"


def _resolve_boundary_condition(
    condition: Mapping[str, str], substitutions: Mapping[str, object]
) -> Dict[str, str]:
    return {
        key: value.format_map(substitutions) if isinstance(value, str) else value
        for key, value in condition.items()
    }


def _render_patch(name: str, condition: Mapping[str, str]) -> str:
    lines = [f"    {name}", "    {", f"        type    {condition['type']};"]
    if "value" in condition:
        lines.append(f"        value   {condition['value']};")
    lines.append("    }")
    return "\n".join(lines)


COUPLING_PROPERTIES_TEMPLATE = """\
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
    class       dictionary;
    location    "constant";
    object      couplingProperties;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Elmer is re-solved when U, T or p has changed by more than these relative
// tolerances since the last update, or after maxStepsBetweenUpdates steps
// (0 = no limit). Zero tolerances update every time step.
velocityTolerance       {velocity_tolerance:g};
temperatureTolerance    {temperature_tolerance:g};
pressureTolerance       {pressure_tolerance:g};
maxStepsBetweenUpdates  {max_steps_between_updates:d};

// Absolute pressure at OpenFOAM p = 0 [Pa], for the relative pressure change
referencePressure       {reference_pressure:g};

// ************************************************************************* //
"""


def render_coupling_properties(config: CaseConfig) -> str:
    """Render constant/couplingProperties for the OpenFOAM solver."""
    values = asdict(config.coupling)
    values["reference_pressure"] = config.plasma.reference_pressure
    return COUPLING_PROPERTIES_TEMPLATE.format_map(values)


class OpenFoamCaseRenderer:
    """Render all OpenFOAM initial and boundary fields."""

    def render(self, config: CaseConfig) -> Dict[str, str]:
        substitutions = {
            "inlet_velocity": _format_vector(config.physics.inlet_velocity),
            "inlet_temperature": f"{config.physics.inlet_temperature:g}",
            "B_field": _format_vector(config.physics.B_field),
        }
        result: Dict[str, str] = {}
        for definition in FIELD_DEFINITIONS:
            result[str(definition["name"])] = self._render_field(
                definition, electrode_patch_names(len(config.electrodes.pairs)), substitutions
            )
        return result

    @staticmethod
    def _render_field(
        definition: Mapping[str, object],
        electrodes: Iterable[str],
        substitutions: Mapping[str, object],
    ) -> str:
        patches: List[str] = []
        for name, condition_name in (
            (INLET, "inlet"),
            (OUTLET, "outlet"),
            (INSULATOR, "wall"),
        ):
            condition = definition[condition_name]
            assert isinstance(condition, Mapping)
            patches.append(
                _render_patch(
                    name, _resolve_boundary_condition(condition, substitutions)
                )
            )
        wall = definition["wall"]
        assert isinstance(wall, Mapping)
        for name in electrodes:
            patches.append(_render_patch(name, _resolve_boundary_condition(wall, substitutions)))
        patches.append(
            _render_patch(
                "defaultFaces", _resolve_boundary_condition(wall, substitutions)
            )
        )

        comment = definition.get("comment")
        return OPENFOAM_HEADER.format(
            field_class=definition["class"],
            field_name=definition["name"],
            comment=f"{comment}\n" if comment else "",
            dimensions=definition["dimensions"],
            internal_field=str(definition["internal_field"]).format_map(substitutions),
            patches="\n\n".join(patches) + "\n",
        )
