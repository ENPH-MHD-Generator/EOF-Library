"""Typed configuration model and validation for linear Hall cases."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import isfinite
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import yaml


class ConfigError(ValueError):
    """Raised when a case configuration is incomplete or physically invalid."""


def _mapping(value: Any, location: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ConfigError(f"{location} must be a mapping")
    return value


def _known_keys(data: Mapping[str, Any], allowed: Iterable[str], location: str) -> None:
    unknown = sorted(set(data) - set(allowed))
    if unknown:
        raise ConfigError(f"Unknown key(s) in {location}: {', '.join(unknown)}")


def _number(value: Any, location: str, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigError(f"{location} must be a number")
    result = float(value)
    if not isfinite(result):
        raise ConfigError(f"{location} must be finite")
    if positive and result <= 0:
        raise ConfigError(f"{location} must be greater than zero")
    return result


def _integer(value: Any, location: str, *, positive: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigError(f"{location} must be an integer")
    if positive and value <= 0:
        raise ConfigError(f"{location} must be greater than zero")
    return value


def _nonnegative_number(value: Any, location: str) -> float:
    result = _number(value, location)
    if result < 0:
        raise ConfigError(f"{location} cannot be negative")
    return result


def _vector(value: Any, location: str) -> Tuple[float, float, float]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != 3:
        raise ConfigError(f"{location} must contain exactly three numbers")
    components = tuple(
        _number(component, f"{location}[{index}]")
        for index, component in enumerate(value)
    )
    return components  # type: ignore[return-value]


@dataclass(frozen=True)
class ChannelConfig:
    length: float
    height: float
    width: float
    wall_thickness: float

    @classmethod
    def from_mapping(cls, raw: Any) -> "ChannelConfig":
        data = _mapping(raw, "channel")
        _known_keys(data, ("length", "height", "width", "wall_thickness"), "channel")
        missing = [
            key
            for key in ("length", "height", "width", "wall_thickness")
            if key not in data
        ]
        if missing:
            raise ConfigError(f"Missing channel key(s): {', '.join(missing)}")
        return cls(
            length=_number(data["length"], "channel.length", positive=True),
            height=_number(data["height"], "channel.height", positive=True),
            width=_number(data["width"], "channel.width", positive=True),
            wall_thickness=_number(
                data["wall_thickness"], "channel.wall_thickness", positive=True
            ),
        )


MESH_TYPES = ("structured", "tetrahedral")
STRUCTURED_MESH_KEYS = ("cell_size", "wall_cell_size", "electrode_edge_cell_size", "growth_rate")
TETRAHEDRAL_MESH_KEYS = ("size_min", "size_max", "size_factor")


@dataclass(frozen=True)
class MeshConfig:
    """Mesh settings.

    ``tetrahedral`` (default): the unstructured Gmsh mesh controlled by the
    size_* keys. ``structured``: graded hexahedra, finest at the walls (cold
    thermal boundary layer) and streamwise at the electrode edges. With the
    local electron energy balance (plasma.electron_energy_transport false) the
    two-temperature model does not converge on structured meshes: their wall
    cells are far smaller than the ~1 mm electron energy relaxation length,
    where the purely local Te balance runs away at electrode edges.
    """

    type: str = "tetrahedral"
    # structured [m]
    cell_size: float = 0.0025  # core cell size
    wall_cell_size: float = 0.00025  # first cell at every wall
    electrode_edge_cell_size: float = 0.001  # streamwise size at electrode edges
    growth_rate: float = 1.2  # maximum size ratio of neighbouring cells
    # tetrahedral
    size_min: Optional[float] = None
    size_max: Optional[float] = None
    size_factor: float = 1.0

    @classmethod
    def from_mapping(cls, raw: Any) -> "MeshConfig":
        data = _mapping(raw or {}, "mesh")
        _known_keys(data, ("type",) + STRUCTURED_MESH_KEYS + TETRAHEDRAL_MESH_KEYS, "mesh")
        mesh_type = data.get("type", cls.type)
        if mesh_type not in MESH_TYPES:
            raise ConfigError(f"mesh.type must be one of: {', '.join(MESH_TYPES)}")
        foreign = STRUCTURED_MESH_KEYS if mesh_type == "tetrahedral" else TETRAHEDRAL_MESH_KEYS
        given = [key for key in foreign if key in data]
        if given:
            raise ConfigError(
                f"mesh key(s) {', '.join(given)} do not apply to type '{mesh_type}'"
                + ("; set 'type: structured' to use them" if mesh_type == "tetrahedral" else "")
            )

        if mesh_type == "structured":
            values = {
                key: _number(data.get(key, getattr(cls, key)), f"mesh.{key}", positive=True)
                for key in STRUCTURED_MESH_KEYS
            }
            if values["growth_rate"] <= 1.0:
                raise ConfigError("mesh.growth_rate must be greater than 1")
            for key in ("wall_cell_size", "electrode_edge_cell_size"):
                if values[key] > values["cell_size"]:
                    raise ConfigError(f"mesh.{key} cannot be larger than mesh.cell_size")
            return cls(type=mesh_type, **values)

        size_min = (
            None
            if data.get("size_min") is None
            else _number(data["size_min"], "mesh.size_min", positive=True)
        )
        size_max = (
            None
            if data.get("size_max") is None
            else _number(data["size_max"], "mesh.size_max", positive=True)
        )
        size_factor = _number(data.get("size_factor", 1.0), "mesh.size_factor", positive=True)
        if size_min is not None and size_max is not None and size_min > size_max:
            raise ConfigError("mesh.size_min cannot be greater than mesh.size_max")
        return cls(
            type=mesh_type, size_min=size_min, size_max=size_max, size_factor=size_factor
        )


@dataclass(frozen=True)
class ElectrodePair:
    x_center: float
    resistance: float


@dataclass(frozen=True)
class ElectrodeConfig:
    length: float
    pairs: Tuple[ElectrodePair, ...]

    @classmethod
    def from_mapping(cls, raw: Any, channel: ChannelConfig) -> "ElectrodeConfig":
        data = _mapping(raw, "electrodes")
        _known_keys(data, ("length", "resistance", "count", "pairs"), "electrodes")

        has_pairs = "pairs" in data
        has_count = "count" in data
        if has_pairs == has_count:
            raise ConfigError("electrodes must define exactly one of 'pairs' or 'count'")

        configured_length = data.get("length")
        if configured_length is not None:
            configured_length = _number(configured_length, "electrodes.length", positive=True)

        if configured_length is None:
            raise ConfigError("electrodes.length is required")

        pairs: List[ElectrodePair] = []
        if has_pairs:
            if "resistance" in data:
                raise ConfigError(
                    "electrodes.resistance is only valid with count; put resistance on each pair"
                )
            raw_pairs = data["pairs"]
            if (
                isinstance(raw_pairs, (str, bytes))
                or not isinstance(raw_pairs, Sequence)
                or not raw_pairs
            ):
                raise ConfigError("electrodes.pairs must be a non-empty list")
            for index, raw_pair in enumerate(raw_pairs, start=1):
                location = f"electrodes.pairs[{index - 1}]"
                pair = _mapping(raw_pair, location)
                _known_keys(pair, ("x_center", "resistance"), location)
                for required in ("x_center", "resistance"):
                    if required not in pair:
                        raise ConfigError(f"{location} is missing '{required}'")
                pairs.append(
                    ElectrodePair(
                        x_center=_number(pair["x_center"], f"{location}.x_center"),
                        resistance=_nonnegative_number(
                            pair["resistance"], f"{location}.resistance"
                        ),
                    )
                )
            length = configured_length
        else:
            count = _integer(data["count"], "electrodes.count", positive=True)
            resistance = _nonnegative_number(
                data.get("resistance", 1.0), "electrodes.resistance"
            )
            spacing = channel.length / (count + 1)
            length = configured_length
            pairs = [
                ElectrodePair(x_center=(index + 1) * spacing, resistance=resistance)
                for index in range(count)
            ]

        cls._validate_geometry(length, pairs, channel)
        return cls(length=length, pairs=tuple(pairs))

    @staticmethod
    def _validate_geometry(
        length: float, pairs: Sequence[ElectrodePair], channel: ChannelConfig
    ) -> None:
        half_length = length / 2.0
        for index, pair in enumerate(pairs, start=1):
            if pair.x_center - half_length < 0 or pair.x_center + half_length > channel.length:
                raise ConfigError(
                    f"Electrode pair {index} extends outside the channel: "
                    f"center={pair.x_center:g}, length={length:g}, "
                    f"channel.length={channel.length:g}"
                )

        ordered = sorted(enumerate(pairs, start=1), key=lambda item: item[1].x_center)
        for (left_index, left), (right_index, right) in zip(ordered, ordered[1:]):
            if right.x_center - left.x_center < length:
                raise ConfigError(
                    f"Electrode pairs {left_index} and {right_index} overlap "
                    f"for shared length {length:g}"
                )


@dataclass(frozen=True)
class PhysicsConfig:
    B_field: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    inlet_velocity: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    inlet_temperature: float = 300.0
    # Fixed wall temperatures [K]; None makes that surface adiabatic. Walls are
    # uncooled and start at room temperature; over a <10 s run a copper
    # electrode surface warms ~10 K and a ceramic one tens to a few hundred K.
    insulator_wall_temperature: Optional[float] = 300.0
    electrode_wall_temperature: Optional[float] = 300.0
    # Absolute static pressure at the outlet [Pa]; the flow is compressible
    outlet_pressure: float = 101325.0

    @classmethod
    def from_mapping(cls, raw: Any) -> "PhysicsConfig":
        data = _mapping(raw or {}, "physics")
        _known_keys(
            data,
            (
                "B_field",
                "inlet_velocity",
                "inlet_temperature",
                "insulator_wall_temperature",
                "electrode_wall_temperature",
                "outlet_pressure",
            ),
            "physics",
        )

        def wall_temperature(key: str) -> Optional[float]:
            value = data.get(key, getattr(cls, key))
            if value is None:
                return None
            return _number(value, f"physics.{key}", positive=True)

        return cls(
            insulator_wall_temperature=wall_temperature("insulator_wall_temperature"),
            electrode_wall_temperature=wall_temperature("electrode_wall_temperature"),
            B_field=_vector(data.get("B_field", (0, 0, 0)), "physics.B_field"),
            inlet_velocity=_vector(
                data.get("inlet_velocity", (0, 0, 0)), "physics.inlet_velocity"
            ),
            inlet_temperature=_number(
                data.get("inlet_temperature", 300),
                "physics.inlet_temperature",
                positive=True,
            ),
            outlet_pressure=_number(
                data.get("outlet_pressure", cls.outlet_pressure),
                "physics.outlet_pressure",
                positive=True,
            ),
        )


@dataclass(frozen=True)
class PlasmaConfig:
    """Alkali-seeded carrier gas; defaults are potassium seed in argon.

    Only the seed ionizes (Saha equation at the electron temperature).
    Conductivity comes from electron-neutral collisions with both the carrier
    gas and neutral seed. With two_temperature, Te is raised above the gas
    temperature by Joule heating of the electrons; otherwise Te = Tgas.

    electron_energy_transport selects how Te is found: true solves the electron
    energy equation (Joule heating, collisional loss to the gas, conduction,
    convection with the gas and the current, and transport of ionization
    energy); false uses the local balance of heating and collisional loss
    (Kerrebrock), which ignores all transport.
    """

    seed_mole_fraction: float = 0.01  # seed atoms per heavy particle
    seed_ionization_energy: float = 4.3407  # eV (K)
    seed_gi_over_gn: float = 0.5  # g(K+) / g(K) = 1 / 2
    seed_cross_section: float = 4.0e-18  # m^2, electron-seed momentum transfer
    carrier_cross_section: float = 1.0e-19  # m^2, electron-carrier momentum transfer
    sigma_min: float = 1.0e-2  # S/m, conductivity floor for matrix conditioning
    sigma_max: float = 1.0e6  # S/m
    two_temperature: bool = True
    electron_energy_transport: bool = True  # electron energy PDE (false: local balance)
    carrier_molar_mass: float = 39.948  # g/mol (Ar)
    seed_molar_mass: float = 39.098  # g/mol (K)
    energy_loss_factor: float = 1.0  # delta; 1 for elastic losses in monatomic gas
    electron_temperature_max: float = 20000.0  # K, cap on the energy-balance solution
    electron_temperature_relaxation: float = 0.5  # under-relaxation of Te per iteration

    @classmethod
    def from_mapping(cls, raw: Any) -> "PlasmaConfig":
        data = _mapping(raw or {}, "plasma")
        defaults = cls()
        if "reference_pressure" in data:
            raise ConfigError(
                "plasma.reference_pressure was removed: the flow is compressible and "
                "OpenFOAM's pressure is absolute; set physics.outlet_pressure instead"
            )
        _known_keys(data, asdict(defaults), "plasma")

        def value(key: str, *, positive: bool = True) -> float:
            location = f"plasma.{key}"
            raw_value = data.get(key, getattr(defaults, key))
            # YAML 1.1 reads exponents without a sign (1.0e6) as strings
            if isinstance(raw_value, str):
                try:
                    raw_value = float(raw_value)
                except ValueError:
                    pass
            if positive:
                return _number(raw_value, location, positive=True)
            return _nonnegative_number(raw_value, location)

        switches = {}
        for key in ("two_temperature", "electron_energy_transport"):
            switches[key] = data.get(key, getattr(defaults, key))
            if not isinstance(switches[key], bool):
                raise ConfigError(f"plasma.{key} must be true or false")

        result = cls(
            seed_mole_fraction=value("seed_mole_fraction"),
            seed_ionization_energy=value("seed_ionization_energy"),
            seed_gi_over_gn=value("seed_gi_over_gn"),
            seed_cross_section=value("seed_cross_section", positive=False),
            carrier_cross_section=value("carrier_cross_section", positive=False),
            sigma_min=value("sigma_min"),
            sigma_max=value("sigma_max"),
            two_temperature=switches["two_temperature"],
            electron_energy_transport=switches["electron_energy_transport"],
            carrier_molar_mass=value("carrier_molar_mass"),
            seed_molar_mass=value("seed_molar_mass"),
            energy_loss_factor=value("energy_loss_factor"),
            electron_temperature_max=value("electron_temperature_max"),
            electron_temperature_relaxation=value("electron_temperature_relaxation"),
        )
        if result.seed_mole_fraction >= 1:
            raise ConfigError("plasma.seed_mole_fraction must be less than 1")
        if result.seed_cross_section + result.carrier_cross_section <= 0:
            raise ConfigError(
                "plasma.seed_cross_section and plasma.carrier_cross_section cannot both be zero"
            )
        if result.sigma_min >= result.sigma_max:
            raise ConfigError("plasma.sigma_min must be less than plasma.sigma_max")
        if result.electron_temperature_relaxation > 1:
            raise ConfigError("plasma.electron_temperature_relaxation must be in (0, 1]")
        return result


@dataclass(frozen=True)
class CouplingConfig:
    """When OpenFOAM re-solves the electrical problem in Elmer.

    The electrical problem is quasi-static: the current depends only on the
    instantaneous velocity, temperature and pressure. Elmer is re-solved once
    any of them has changed by more than its relative tolerance since the last
    update, or after max_steps_between_updates steps. All zero (the default)
    updates every time step.
    """

    velocity_tolerance: float = 0.0  # max |U - U_sent| / max |U_sent|
    temperature_tolerance: float = 0.0  # Joule-power-weighted RMS of |T - T_sent| / T_sent
    pressure_tolerance: float = 0.0  # max |p - p_sent| / absolute pressure
    max_steps_between_updates: int = 0  # 0 = no limit

    @classmethod
    def from_mapping(cls, raw: Any) -> "CouplingConfig":
        data = _mapping(raw or {}, "coupling")
        defaults = cls()
        _known_keys(data, asdict(defaults), "coupling")
        return cls(
            velocity_tolerance=_nonnegative_number(
                data.get("velocity_tolerance", defaults.velocity_tolerance),
                "coupling.velocity_tolerance",
            ),
            temperature_tolerance=_nonnegative_number(
                data.get("temperature_tolerance", defaults.temperature_tolerance),
                "coupling.temperature_tolerance",
            ),
            pressure_tolerance=_nonnegative_number(
                data.get("pressure_tolerance", defaults.pressure_tolerance),
                "coupling.pressure_tolerance",
            ),
            max_steps_between_updates=_integer(
                data.get("max_steps_between_updates", defaults.max_steps_between_updates),
                "coupling.max_steps_between_updates",
            ),
        )


LINEAR_SOLVERS = ("auto", "iterative", "mumps")


@dataclass(frozen=True)
class NumericsConfig:
    """Numerical method choices.

    ``linear_solver`` for Elmer's potential equation: ``iterative`` (ILU-
    preconditioned GCR), ``mumps`` (parallel sparse direct), or ``auto``, which
    uses iterative on tetrahedral meshes (faster there) and MUMPS on structured
    meshes, where ILU converges slowly on the thin wall cells.
    """

    linear_solver: str = "auto"

    @classmethod
    def from_mapping(cls, raw: Any) -> "NumericsConfig":
        data = _mapping(raw or {}, "numerics")
        _known_keys(data, ("linear_solver",), "numerics")
        solver = data.get("linear_solver", cls.linear_solver)
        if solver not in LINEAR_SOLVERS:
            raise ConfigError(f"numerics.linear_solver must be one of: {', '.join(LINEAR_SOLVERS)}")
        return cls(linear_solver=solver)


@dataclass(frozen=True)
class CaseConfig:
    channel: ChannelConfig
    mesh: MeshConfig
    electrodes: ElectrodeConfig
    physics: PhysicsConfig
    plasma: PlasmaConfig = PlasmaConfig()
    coupling: CouplingConfig = CouplingConfig()
    numerics: NumericsConfig = NumericsConfig()
    schema_version: int = 1

    @classmethod
    def from_mapping(cls, raw: Any) -> "CaseConfig":
        data = _mapping(raw, "configuration")
        _known_keys(
            data,
            ("schema_version", "channel", "mesh", "electrodes", "physics", "plasma", "coupling", "numerics"),
            "configuration",
        )
        schema_version = _integer(data.get("schema_version", 1), "schema_version", positive=True)
        if schema_version != 1:
            raise ConfigError(f"Unsupported schema_version {schema_version}; expected 1")
        if "channel" not in data or "electrodes" not in data:
            raise ConfigError("configuration requires 'channel' and 'electrodes' sections")
        channel = ChannelConfig.from_mapping(data["channel"])
        return cls(
            channel=channel,
            mesh=MeshConfig.from_mapping(data.get("mesh", {})),
            electrodes=ElectrodeConfig.from_mapping(data["electrodes"], channel),
            physics=PhysicsConfig.from_mapping(data.get("physics", {})),
            plasma=PlasmaConfig.from_mapping(data.get("plasma", {})),
            coupling=CouplingConfig.from_mapping(data.get("coupling", {})),
            numerics=NumericsConfig.from_mapping(data.get("numerics", {})),
            schema_version=schema_version,
        )

    def to_dict(self) -> Dict[str, Any]:
        """Return the canonical, fully resolved public representation."""
        return {
            "schema_version": self.schema_version,
            "channel": asdict(self.channel),
            "mesh": {
                key: value
                for key, value in asdict(self.mesh).items()
                if key == "type"
                or key
                in (
                    STRUCTURED_MESH_KEYS
                    if self.mesh.type == "structured"
                    else TETRAHEDRAL_MESH_KEYS
                )
            },
            "electrodes": {
                "length": self.electrodes.length,
                "pairs": [asdict(pair) for pair in self.electrodes.pairs],
            },
            "physics": {
                "B_field": list(self.physics.B_field),
                "inlet_velocity": list(self.physics.inlet_velocity),
                "inlet_temperature": self.physics.inlet_temperature,
                "insulator_wall_temperature": self.physics.insulator_wall_temperature,
                "electrode_wall_temperature": self.physics.electrode_wall_temperature,
                "outlet_pressure": self.physics.outlet_pressure,
            },
            "plasma": asdict(self.plasma),
            "coupling": asdict(self.coupling),
            "numerics": asdict(self.numerics),
        }


def _load_yaml_mapping(path: Path) -> Dict[str, Any]:
    try:
        with path.open(encoding="utf-8") as stream:
            raw = yaml.safe_load(stream)
    except OSError as exc:
        raise ConfigError(f"Cannot read configuration '{path}': {exc}") from exc
    except yaml.YAMLError as exc:
        raise ConfigError(f"Invalid YAML in '{path}': {exc}") from exc
    if not isinstance(raw, dict):
        raise ConfigError(f"Configuration '{path}' must contain a YAML mapping")
    return raw


def load_case_config(path: Path) -> CaseConfig:
    return CaseConfig.from_mapping(_load_yaml_mapping(path))
