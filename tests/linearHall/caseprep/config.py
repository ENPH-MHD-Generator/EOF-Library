"""Compatibility imports for the configuration API.

New code should import runtime models from :mod:`caseprep.models` and the YAML
loader from :mod:`caseprep.loader`.
"""

from .loader import ConfigError, load_case_config
from .models import (
    CaseConfig,
    CaseConfigV1,
    ChannelConfig,
    ElectrodeConfig,
    ElectrodePair,
    EvenlySpacedElectrodeConfig,
    ExplicitElectrodeConfig,
    MeshConfig,
    PhysicsConfig,
)

__all__ = [
    "CaseConfig",
    "CaseConfigV1",
    "ChannelConfig",
    "ConfigError",
    "ElectrodeConfig",
    "ElectrodePair",
    "EvenlySpacedElectrodeConfig",
    "ExplicitElectrodeConfig",
    "MeshConfig",
    "PhysicsConfig",
    "load_case_config",
]
