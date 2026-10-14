"""Canonical Pydantic models for MHD experiments."""

from __future__ import annotations

from math import isfinite
from pathlib import PurePosixPath
from typing import Annotated, Any, Dict, List, Literal, Mapping, Optional, Tuple, Union

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

    insulator_wall_temperature: Optional[PositiveNumber] = 300.0
    electrode_wall_temperature: Optional[PositiveNumber] = 300.0
    outlet_pressure: PositiveNumber = 101325.0


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


class TetrahedralMeshConfig(SchemaModel):
    """Uniform target length in metres; the sole tetrahedral sizing policy."""

    target_element_size: PositiveNumber = Field(
        description=(
            "Uniform target characteristic length for the shared Gmsh "
            "tetrahedral mesh."
        ),
        json_schema_extra={"x-unit": "m"},
    )

    type: Literal["tetrahedral"] = "tetrahedral"


class StructuredMeshConfig(SchemaModel):
    """Graded hexahedra with separate wall, edge, and core resolution."""

    type: Literal["structured"]
    cell_size: PositiveNumber = 0.003
    wall_cell_size: PositiveNumber = 0.0005
    electrode_edge_cell_size: PositiveNumber = 0.0005
    growth_rate: Annotated[FiniteNumber, Field(gt=1)] = 1.2
    electrode_wall_cell_size: Optional[PositiveNumber] = 1.5e-5
    side_wall_cell_size: Optional[PositiveNumber] = None
    streamwise_cell_size: Optional[PositiveNumber] = 0.005

    @model_validator(mode="after")
    def validate_grading(self) -> "StructuredMeshConfig":
        for name in ("wall_cell_size", "electrode_wall_cell_size", "side_wall_cell_size"):
            size = getattr(self, name)
            if size is None:
                size = self.wall_cell_size
            if size > self.cell_size:
                raise ValueError(f"mesh.{name} cannot be larger than mesh.cell_size")
        streamwise = self.streamwise_cell_size or self.cell_size
        if self.electrode_edge_cell_size > streamwise:
            raise ValueError("mesh.electrode_edge_cell_size cannot be larger than mesh.streamwise_cell_size")
        return self


MeshConfig = Annotated[
    Union[TetrahedralMeshConfig, StructuredMeshConfig], Field(discriminator="type")
]


class PlasmaConfig(SchemaModel):
    """Seeded argon plasma, with collision data supplied by plasma_collisions."""

    seed_mole_fraction: Annotated[PositiveNumber, Field(lt=1)] = 0.01
    electron_transport_model: Literal["drifting", "lorentz"] = "drifting"
    potassium_elastic_scale: PositiveNumber = 1.0
    excited_state_temperature: Literal["electron", "gas"] = "electron"
    ion_reduced_mobility: NonNegativeNumber = 2.43e-4
    electron_wall_heat_transmission: NonNegativeNumber = 0.0
    sigma_min: PositiveNumber = 1.0e-2
    sigma_max: PositiveNumber = 1.0e6
    two_temperature: Annotated[bool, Field(strict=True)] = True
    electron_energy_transport: Annotated[bool, Field(strict=True)] = True
    electron_temperature_max: PositiveNumber = 20000.0
    electron_temperature_relaxation: Annotated[PositiveNumber, Field(le=1)] = 0.5

    @model_validator(mode="after")
    def validate_conductivity_bounds(self) -> "PlasmaConfig":
        if self.sigma_min >= self.sigma_max:
            raise ValueError("plasma.sigma_min must be less than plasma.sigma_max")
        return self


class CouplingConfig(SchemaModel):
    """Relative state-change thresholds; all zero updates every time step."""

    velocity_tolerance: NonNegativeNumber = 0.0
    temperature_tolerance: NonNegativeNumber = 0.0
    pressure_tolerance: NonNegativeNumber = 0.0
    max_steps_between_updates: Annotated[int, Field(strict=True, ge=0)] = 0


class NumericsConfig(SchemaModel):
    linear_solver: Literal["auto", "iterative", "mumps"] = "auto"


class CaseConfig(SchemaModel):
    """Canonical configuration for the compressible Saha MHD pipeline."""

    model_config = ConfigDict(title="MHD Experiment Configuration")

    channel: ChannelConfig
    mesh: MeshConfig
    electrodes: ElectrodeConfig
    physics: PhysicsConfig = Field(default_factory=PhysicsConfig)
    plasma: PlasmaConfig = Field(default_factory=PlasmaConfig)
    coupling: CouplingConfig = Field(default_factory=CouplingConfig)
    numerics: NumericsConfig = Field(default_factory=NumericsConfig)
    openfoam: OpenFoamConfig

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
    def validate_electrode_geometry(self) -> "CaseConfig":
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

    @model_validator(mode="after")
    def validate_solver_inputs(self) -> "CaseConfig":
        invalid_files = set(self.openfoam.files) & {"constant/transportProperties", "constant/g", "0/p_rgh"}
        if invalid_files:
            raise ValueError("The compressible Saha solver does not use: " + ", ".join(sorted(invalid_files)))
        if "constant/couplingProperties" in self.openfoam.files:
            raise ValueError("Set coupling properties through the canonical coupling section")
        return self

    def to_dict(self) -> Dict[str, Any]:
        """Serialize defaults and explicit electrode pairs for a prepared case."""
        result = self.model_dump(mode="json", by_alias=True)
        result["electrodes"] = {
            "length": self.electrodes.length,
            "pairs": [pair.model_dump(mode="json") for pair in self.electrode_pairs],
        }
        return result
