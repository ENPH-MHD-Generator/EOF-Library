"""Electron collision data for Elmer, from `plasma_collisions`.

Elmer evaluates the plasma state at every node on every nonlinear iteration, so
it cannot call Python. Case preparation tabulates everything that depends on
the electron temperature alone (and, for the Lorentz model, on two more reduced
variables) into a plain-text file that MHDSolve reads once and interpolates.
What depends on the local state (electron density, Coulomb logarithm, gas and
electron temperature, magnetic field) is assembled in Fortran.

All coefficients come from a Maxwellian EEDF at `theta_e` for argon seeded with
potassium, with the cross sections `plasma_collisions` provides. They are built
for a mixture at a negligible gas temperature with superelastic collisions off,
which leaves the pure electron-temperature integrals; MHDSolve applies the gas
temperature (recoil factor `1 - Tg/Te`, excited-state populations) itself.

File format, read by `ReadCollisionData` in MHDUtils.F90: lines starting with
``!`` are comments; then blocks of ``<name> <count>`` followed by `count`
whitespace-separated values; ``end`` closes the file. Three-dimensional blocks
are written in Fortran (column-major) order: theta varies fastest.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np

from .config import CaseConfig

KELVIN_PER_EV = 1.602_176_634e-19 / 1.380_649e-23
ELECTRON_CHARGE_OVER_MASS = 1.602_176_634e-19 / 9.109_383_7e-31

TABLE_FILE_NAME = "electron_collisions.dat"
METADATA_FILE_NAME = "electron_collisions.json"
FORMAT_VERSION = 1

# Reduced-variable grids for the Lorentz model. The Coulomb parameter is
# n_e ln(Lambda) / N (the electron-ion momentum cross section is this times a
# fixed function of energy); the magnetic parameter is omega_ce / N. Below the
# lower ends both effects are negligible and MHDSolve clamps to them.
COULOMB_PARAMETER_RANGE = (1e-10, 1.0)
MAGNETIC_PARAMETER_RANGE = (1e-19, 1e-8)  # m^3/s: 1.5 T at 3e24 m^-3 is ~9e-14


@dataclass(frozen=True)
class CollisionTables:
    """Everything written to the Elmer table file."""

    blocks: Dict[str, np.ndarray]
    metadata: Dict[str, object]
    saha: Dict[str, float]
    """Seed ionization energy [eV] and statistical weight ratio g_ion/g_neutral."""


def _mixture(config: CaseConfig, *, gas_temperature: float, superelastic: bool):
    import plasma_collisions as pc

    return pc.argon_potassium(
        seed_fraction=config.plasma.seed_mole_fraction,
        gas_temperature=gas_temperature,
        excitation_temperature=None if superelastic else 0.0,
        potassium_elastic_scale=config.plasma.potassium_elastic_scale,
    )


def theta_grid(config: CaseConfig, points: int = 90) -> np.ndarray:
    """Electron temperatures to tabulate, eV: from below the coldest wall to past Te max."""
    upper = max(3.0, 1.1 * config.plasma.electron_temperature_max / KELVIN_PER_EV)
    return np.geomspace(0.015, upper, points)


def build_tables(config: CaseConfig) -> CollisionTables:
    """Tabulate the electron coefficients MHDSolve needs for this case."""
    import plasma_collisions as pc

    with np.errstate(all="ignore"):  # the library's quadrature over/underflows harmlessly
        return _build_tables(config, pc)


def _build_tables(config: CaseConfig, pc) -> CollisionTables:
    # Pure electron-temperature integrals: negligible gas temperature (recoil
    # factor 1), no superelastic collisions. See the module docstring.
    cold = _mixture(config, gas_temperature=1e-3, superelastic=False)
    rates = pc.maxwellian(cold)
    theta = theta_grid(config)
    potassium = cold["K"]
    blocks: Dict[str, np.ndarray] = {"theta": theta}

    # Drifting Maxwellian: electron-neutral momentum-transfer frequency / N. The
    # NRL electron-ion frequency adds to it exactly, and is evaluated in Fortran.
    blocks["momentum_frequency_n"] = rates.momentum_frequency_n(theta)
    # Recoil energy loss to neutrals / N without the (1 - Tg/Te) factor, eV m^3/s.
    blocks["elastic_recoil_n"] = rates.elastic_power_loss_n(theta)
    # Recoil loss to K+ per unit Coulomb parameter a = n_e ln(Lambda)/N, eV m^3/s.
    # transport() returns e n_e N a G_ei for a gas this cold; divide it back out.
    n_ref, a_ref = 1e24, 1e-6
    probe = rates.transport(
        theta, n_ref, a_ref * n_ref, 0.0, model="drifting", ln_lambda=1.0
    )
    watts = 1.602_176_634e-19 * (a_ref * n_ref) * n_ref
    blocks["ei_recoil_per_coulomb"] = np.asarray(probe.ei_power_loss) / watts / a_ref
    # K 4s-4p excitation and its superelastic inverse (per excited atom), m^3/s.
    blocks["k_excitation_K"] = rates.rate_coefficient(theta, "K", "excitation")
    blocks["k_superelastic_K"] = rates.superelastic_rate_coefficient(theta, "K")
    excitation = potassium.processes(pc.ProcessType.EXCITATION)[0]
    blocks["excitation_energy_K"] = np.array([float(excitation.threshold)])
    blocks["excitation_weight_ratio_K"] = np.array([float(excitation.weight_ratio or 1.0)])

    model = config.plasma.electron_transport_model
    if model == "lorentz":
        blocks.update(_lorentz_blocks(rates, theta))

    ground = potassium.ground_weight
    ion_partition = sum(g for energy, g in potassium.ion_levels if energy == 0.0)
    saha = {
        "ionization_energy": float(potassium.ionization_energy),
        "weight_ratio": float(ion_partition / ground),
    }
    metadata: Dict[str, object] = {
        "generator": "caseprep.collisions via plasma_collisions",
        "plasma_collisions_version": _package_version(),
        "eedf": "maxwellian",
        "electron_transport_model": model,
        "potassium_elastic_scale": config.plasma.potassium_elastic_scale,
        "seed_mole_fraction": config.plasma.seed_mole_fraction,
        "sources": cold.metadata()["sources"],
    }
    return CollisionTables(blocks=blocks, metadata=metadata, saha=saha)


def _lorentz_blocks(rates, theta: np.ndarray) -> Dict[str, np.ndarray]:
    """Mobility x N, Pedersen and Hall, on (theta, Coulomb parameter, magnetic parameter)."""
    # Trilinear interpolation in the logarithms against the library: Pedersen
    # and Hall mobilities within 2% (95th percentile) and 3% (max) over
    # theta 0.02-1.5 eV, x_e 1e-8-1e-2, B 0.05-5 T, N 1e23-3e25 m^-3.
    coulomb = np.geomspace(*COULOMB_PARAMETER_RANGE, 41)
    magnetic = np.geomspace(*MAGNETIC_PARAMETER_RANGE, 61)
    n_ref = 1e24
    # Fortran order: index (i_theta, i_coulomb, i_magnetic), theta fastest.
    t, a, w = np.meshgrid(theta, coulomb, magnetic, indexing="ij")
    result = rates.transport(
        t.ravel(order="F"),
        n_ref,
        a.ravel(order="F") * n_ref,
        w.ravel(order="F") * n_ref / ELECTRON_CHARGE_OVER_MASS,
        model="lorentz",
        ln_lambda=1.0,
    )
    return {
        "coulomb_parameter": coulomb,
        "magnetic_parameter": magnetic,
        "mobility_n": np.asarray(result.mobility) * n_ref,
        "pedersen_mobility_n": np.asarray(result.pedersen_mobility) * n_ref,
        "hall_mobility_n": np.asarray(result.hall_mobility) * n_ref,
    }


def _package_version() -> str:
    try:
        from importlib.metadata import version

        return version("plasma-collisions")
    except Exception:  # pragma: no cover - metadata is informational
        return "unknown"


def render_table_file(tables: CollisionTables) -> str:
    lines: List[str] = [
        "! Electron collision data for MHDSolve, generated by caseprep from plasma_collisions.",
        "! Do not edit; regenerate with mhd prepare. Provenance: electron_collisions.json",
        "format 1",
        str(FORMAT_VERSION),
    ]
    for name, values in tables.blocks.items():
        array = np.asarray(values, dtype=np.float64).ravel()
        lines.append(f"{name} {array.size}")
        lines.extend(_rows(array))
    lines.append("end 0")
    return "\n".join(lines) + "\n"


def _rows(values: Sequence[float], per_line: int = 6) -> List[str]:
    return [
        " ".join(f"{v:.7e}" for v in values[i : i + per_line])
        for i in range(0, len(values), per_line)
    ]


def write_collision_tables(config: CaseConfig, directory: Path) -> CollisionTables:
    """Write the Elmer table file into `directory`."""
    tables = build_tables(config)
    directory.mkdir(parents=True, exist_ok=True)
    (directory / TABLE_FILE_NAME).write_text(render_table_file(tables), encoding="utf-8")
    (directory / METADATA_FILE_NAME).write_text(
        json.dumps({**tables.metadata, "saha": tables.saha}, indent=2) + "\n", encoding="utf-8"
    )
    return tables


def langevin_reduced_mobility(
    polarizability_angstrom3: float = 1.6411, ion_mass_amu: float = 39.0983, gas_mass_amu: float = 39.948
) -> float:
    """Langevin (polarization) reduced mobility of an ion in a gas, m^2/(V s) at N0 = 2.6868e25 m^-3.

    `K0 = 13.853 / sqrt(alpha mu)` cm^2/(V s), alpha in A^3 and the reduced mass
    mu in amu. The defaults are K+ in argon.
    """
    reduced = ion_mass_amu * gas_mass_amu / (ion_mass_amu + gas_mass_amu)
    return 13.853 / np.sqrt(polarizability_angstrom3 * reduced) * 1e-4
