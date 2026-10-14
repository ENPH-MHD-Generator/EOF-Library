"""Versioned Pydantic models for the public MHD experiment schema."""

from __future__ import annotations

from math import isfinite
from pathlib import PurePosixPath
from typing import Annotated, Any, Dict, List, Literal, Mapping, Tuple, Union

from pydantic import (
    BaseModel,
    BeforeValidator,
    ConfigDict,
    Field,
    field_validator,
    model_validator,
)
from typing_extensions import TypeAliasType


def _finite_number(value: Any) -> float:
    """Accept YAML integers and floats, but reject booleans and strings."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("must be a number")
    try:
        result = float(value)
    except OverflowError as exc:
        raise ValueError("must be finite") from exc
    if not isfinite(result):
        raise ValueError("must be finite")
    return result


FiniteNumber = Annotated[float, BeforeValidator(_finite_number)]
PositiveNumber = Annotated[FiniteNumber, Field(gt=0)]
NonNegativeNumber = Annotated[FiniteNumber, Field(ge=0)]
PositiveInteger = Annotated[int, Field(strict=True, gt=0)]
Vector3 = Tuple[FiniteNumber, FiniteNumber, FiniteNumber]
FoamNumber = Annotated[float, Field(allow_inf_nan=False)]
OpenFoamListValue = TypeAliasType(
    "OpenFoamListValue",
    Union[
        bool, int, FoamNumber, str,
        List["OpenFoamListValue"], Dict[str, "OpenFoamValue"],
    ],
)
OpenFoamValue = TypeAliasType(
    "OpenFoamValue",
    Union[
        None,
        bool,
        int,
        FoamNumber,
        str,
        List[OpenFoamListValue],
        Dict[str, "OpenFoamValue"],
    ],
)


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


class OpenFoamFileConfig(SchemaModel):
    """One data-driven OpenFOAM dictionary-like file."""

    model_config = ConfigDict(
        extra="forbid",
        frozen=True,
        validate_default=True,
        populate_by_name=True,
    )

    foam_class: str = Field(
        default="dictionary",
        alias="class",
        description="OpenFOAM FoamFile class.",
    )
    entries: Dict[str, OpenFoamValue] = Field(
        description=(
            "Ordered OpenFOAM entries. Strings contain OpenFOAM tokens, mappings "
            "become sub-dictionaries, arrays become parenthesized lists, booleans "
            "become yes/no, and null emits a directive without a value."
        )
    )

    @field_validator("foam_class")
    @classmethod
    def validate_foam_class(cls, value: str) -> str:
        if not value or any(character.isspace() for character in value):
            raise ValueError("must be a non-empty OpenFOAM class token")
        return value


class OpenFoamConfig(SchemaModel):
    """OpenFOAM files rendered from the central experiment configuration."""

    files: Dict[str, OpenFoamFileConfig] = Field(min_length=1)

    @field_validator("files")
    @classmethod
    def validate_file_paths(
        cls, value: Dict[str, OpenFoamFileConfig]
    ) -> Dict[str, OpenFoamFileConfig]:
        for raw_path in value:
            path = PurePosixPath(raw_path)
            if (
                path.is_absolute()
                or path.as_posix() != raw_path
                or len(path.parts) != 2
                or path.parts[0] not in {"0", "constant", "system"}
                or any(part in {"", ".", ".."} for part in path.parts)
                or any(character.isspace() or character in '\\";{}()' for character in raw_path)
            ):
                raise ValueError(
                    f"'{raw_path}' must be a file directly beneath 0/, constant/, "
                    "or system/"
                )
        return value


class CaseConfigBase(SchemaModel):
    """Shared experiment fields and cross-section physical constraints."""

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
    def validate_electrode_geometry(self) -> "CaseConfigBase":
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


class CaseConfigV1(CaseConfigBase):
    """Legacy schema-version-1 experiment configuration."""

    model_config = ConfigDict(title="MHD Experiment Configuration v1")

    schema_version: Literal[1]


class CaseConfigV2(CaseConfigBase):
    """Experiment with data-driven OpenFOAM dictionary rendering."""

    model_config = ConfigDict(title="MHD Experiment Configuration v2")

    schema_version: Literal[2]
    openfoam: OpenFoamConfig

    def to_dict(self) -> Dict[str, Any]:
        result = super().to_dict()
        result["openfoam"] = self.openfoam.model_dump(mode="json", by_alias=True)
        return result


# Runtime code consumes this stable name. Future versions get separate models;
# caseprep.loader retains explicit schema_version dispatch.
CaseConfig = Union[CaseConfigV1, CaseConfigV2]
