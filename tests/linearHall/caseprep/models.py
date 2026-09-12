"""Versioned Pydantic models for the public MHD experiment schema."""

from __future__ import annotations

from math import isfinite
from typing import Annotated, Any, Dict, Literal, Mapping, Tuple, Union

from pydantic import BaseModel, BeforeValidator, ConfigDict, Field, model_validator


def _finite_number(value: Any) -> float:
    """Accept YAML integers and floats, but reject booleans and strings."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("must be a number")
    result = float(value)
    if not isfinite(result):
        raise ValueError("must be finite")
    return result


FiniteNumber = Annotated[float, BeforeValidator(_finite_number)]
PositiveNumber = Annotated[FiniteNumber, Field(gt=0)]
NonNegativeNumber = Annotated[FiniteNumber, Field(ge=0)]
PositiveInteger = Annotated[int, Field(strict=True, gt=0)]
Vector3 = Tuple[FiniteNumber, FiniteNumber, FiniteNumber]


class SchemaModel(BaseModel):
    """Shared behavior for immutable, closed experiment-schema objects."""

    model_config = ConfigDict(
        extra="forbid",
        frozen=True,
        validate_default=True,
    )


class ChannelConfig(SchemaModel):
    length: PositiveNumber = Field(
        description="Flow-channel length.",
        json_schema_extra={"x-unit": "m"},
    )
    height: PositiveNumber = Field(
        description="Flow-channel height.",
        json_schema_extra={"x-unit": "m"},
    )
    width: PositiveNumber = Field(
        description="Flow-channel width.",
        json_schema_extra={"x-unit": "m"},
    )
    wall_thickness: PositiveNumber = Field(
        description="Electrode and insulating shell thickness.",
        json_schema_extra={"x-unit": "m"},
    )


class MeshConfig(SchemaModel):
    target_element_size: PositiveNumber = Field(
        description=(
            "Uniform target characteristic length for the shared Gmsh "
            "tetrahedral mesh."
        ),
        json_schema_extra={"x-unit": "m"},
    )

    @model_validator(mode="before")
    @classmethod
    def reject_legacy_controls(cls, value: Any) -> Any:
        """Give old prototype configurations an actionable migration error."""
        if isinstance(value, Mapping):
            legacy = sorted(
                set(value) & {"size_min", "size_max", "size_factor"}
            )
            if legacy:
                names = ", ".join(f"mesh.{name}" for name in legacy)
                raise ValueError(
                    f"{names} no longer supported; replace the legacy mesh "
                    "controls with mesh.target_element_size in metres"
                )
        return value


class ElectrodePair(SchemaModel):
    x_center: FiniteNumber = Field(
        description="Electrode-pair center along the channel x-axis.",
        json_schema_extra={"x-unit": "m"},
    )
    resistance: NonNegativeNumber = Field(
        description="External resistance across this electrode pair.",
        json_schema_extra={"x-unit": "ohm"},
    )


class ExplicitElectrodeConfig(SchemaModel):
    """Explicitly positioned electrode pairs with individual resistances."""

    length: PositiveNumber = Field(
        description="Shared streamwise length for every electrode.",
        json_schema_extra={"x-unit": "m"},
    )
    pairs: Annotated[Tuple[ElectrodePair, ...], Field(min_length=1)]


class EvenlySpacedElectrodeConfig(SchemaModel):
    """A requested number of identical, evenly spaced electrode pairs."""

    length: PositiveNumber = Field(
        description="Shared streamwise length for every electrode.",
        json_schema_extra={"x-unit": "m"},
    )
    count: PositiveInteger
    resistance: NonNegativeNumber = Field(
        default=1.0,
        description="External resistance applied to every generated pair.",
        json_schema_extra={"x-unit": "ohm"},
    )


ElectrodeConfig = Union[
    ExplicitElectrodeConfig,
    EvenlySpacedElectrodeConfig,
]


class PhysicsConfig(SchemaModel):
    B_field: Vector3 = Field(
        default=(0.0, 0.0, 0.0),
        description="Applied magnetic-field vector.",
        json_schema_extra={"x-unit": "T"},
    )
    inlet_velocity: Vector3 = Field(
        default=(0.0, 0.0, 0.0),
        description="OpenFOAM inlet velocity vector.",
        json_schema_extra={"x-unit": "m/s"},
    )
    inlet_temperature: PositiveNumber = Field(
        default=300.0,
        description="OpenFOAM inlet temperature.",
        json_schema_extra={"x-unit": "K"},
    )


class CaseConfigV1(SchemaModel):
    """Complete schema-version-1 experiment and its physical constraints."""

    model_config = ConfigDict(title="MHD Experiment Configuration v1")

    schema_version: Literal[1]
    channel: ChannelConfig
    mesh: MeshConfig
    electrodes: ElectrodeConfig
    physics: PhysicsConfig = Field(default_factory=PhysicsConfig)

    @property
    def electrode_pairs(self) -> Tuple[ElectrodePair, ...]:
        """Return explicit pairs regardless of the placement syntax in the YAML."""
        if isinstance(self.electrodes, ExplicitElectrodeConfig):
            return self.electrodes.pairs

        spacing = self.channel.length / (self.electrodes.count + 1)
        return tuple(
            ElectrodePair(
                x_center=(index + 1) * spacing,
                resistance=self.electrodes.resistance,
            )
            for index in range(self.electrodes.count)
        )

    @model_validator(mode="after")
    def validate_electrode_geometry(self) -> "CaseConfigV1":
        """Validate constraints that span channel and electrode sections."""
        length = self.electrodes.length
        pairs = self.electrode_pairs
        half_length = length / 2.0

        for index, pair in enumerate(pairs, start=1):
            if (
                pair.x_center - half_length < 0
                or pair.x_center + half_length > self.channel.length
            ):
                raise ValueError(
                    f"electrode pair {index} extends outside the channel: "
                    f"center={pair.x_center:g}, length={length:g}, "
                    f"channel.length={self.channel.length:g}"
                )

        ordered = sorted(
            enumerate(pairs, start=1),
            key=lambda item: item[1].x_center,
        )
        for (left_index, left), (right_index, right) in zip(
            ordered, ordered[1:]
        ):
            if right.x_center - left.x_center < length:
                raise ValueError(
                    f"electrode pairs {left_index} and {right_index} overlap "
                    f"for shared length {length:g}"
                )
        return self

    def to_dict(self) -> Dict[str, Any]:
        """Return the canonical, fully resolved public representation."""
        return {
            "schema_version": self.schema_version,
            "channel": self.channel.model_dump(mode="json"),
            "mesh": self.mesh.model_dump(mode="json"),
            "electrodes": {
                "length": self.electrodes.length,
                "pairs": [
                    pair.model_dump(mode="json") for pair in self.electrode_pairs
                ],
            },
            "physics": self.physics.model_dump(mode="json"),
        }


# Runtime code consumes this stable name. Future versions get separate models;
# caseprep.loader retains explicit schema_version dispatch.
CaseConfig = CaseConfigV1
