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


@dataclass(frozen=True)
class MeshConfig:
    target_element_size: float

    @classmethod
    def from_mapping(cls, raw: Any) -> "MeshConfig":
        data = _mapping(raw, "mesh")
        legacy_keys = sorted(set(data) & {"size_min", "size_max", "size_factor"})
        if legacy_keys:
            legacy_names = ", ".join(f"mesh.{key}" for key in legacy_keys)
            raise ConfigError(
                f"{legacy_names} no longer supported; replace the legacy mesh "
                "controls with mesh.target_element_size in metres"
            )
        _known_keys(data, ("target_element_size",), "mesh")
        if "target_element_size" not in data:
            raise ConfigError("mesh.target_element_size is required")
        return cls(
            target_element_size=_number(
                data["target_element_size"],
                "mesh.target_element_size",
                positive=True,
            )
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

    @classmethod
    def from_mapping(cls, raw: Any) -> "PhysicsConfig":
        data = _mapping(raw or {}, "physics")
        _known_keys(data, ("B_field", "inlet_velocity", "inlet_temperature"), "physics")
        return cls(
            B_field=_vector(data.get("B_field", (0, 0, 0)), "physics.B_field"),
            inlet_velocity=_vector(
                data.get("inlet_velocity", (0, 0, 0)), "physics.inlet_velocity"
            ),
            inlet_temperature=_number(
                data.get("inlet_temperature", 300),
                "physics.inlet_temperature",
                positive=True,
            ),
        )


@dataclass(frozen=True)
class CaseConfig:
    channel: ChannelConfig
    mesh: MeshConfig
    electrodes: ElectrodeConfig
    physics: PhysicsConfig
    schema_version: int = 1

    @classmethod
    def from_mapping(cls, raw: Any) -> "CaseConfig":
        data = _mapping(raw, "configuration")
        _known_keys(
            data,
            ("schema_version", "channel", "mesh", "electrodes", "physics"),
            "configuration",
        )
        schema_version = _integer(data.get("schema_version", 1), "schema_version", positive=True)
        if schema_version != 1:
            raise ConfigError(f"Unsupported schema_version {schema_version}; expected 1")
        missing_sections = [
            section
            for section in ("channel", "mesh", "electrodes")
            if section not in data
        ]
        if missing_sections:
            raise ConfigError(
                "configuration requires section(s): " + ", ".join(missing_sections)
            )
        channel = ChannelConfig.from_mapping(data["channel"])
        return cls(
            channel=channel,
            mesh=MeshConfig.from_mapping(data["mesh"]),
            electrodes=ElectrodeConfig.from_mapping(data["electrodes"], channel),
            physics=PhysicsConfig.from_mapping(data.get("physics", {})),
            schema_version=schema_version,
        )

    def to_dict(self) -> Dict[str, Any]:
        """Return the canonical, fully resolved public representation."""
        return {
            "schema_version": self.schema_version,
            "channel": asdict(self.channel),
            "mesh": asdict(self.mesh),
            "electrodes": {
                "length": self.electrodes.length,
                "pairs": [asdict(pair) for pair in self.electrodes.pairs],
            },
            "physics": {
                "B_field": list(self.physics.B_field),
                "inlet_velocity": list(self.physics.inlet_velocity),
                "inlet_temperature": self.physics.inlet_temperature,
            },
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
