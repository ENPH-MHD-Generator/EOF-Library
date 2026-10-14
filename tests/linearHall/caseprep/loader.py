"""YAML loading and user-facing validation errors for the canonical schema."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Mapping, Tuple, Union

import yaml
from pydantic import ValidationError

from .models import CaseConfig


class ConfigError(ValueError):
    """Raised when an experiment file cannot be loaded or validated."""



class UniqueKeyLoader(yaml.SafeLoader):
    """Reject ambiguous YAML rather than silently using the last value."""

    def construct_mapping(self, node: yaml.MappingNode, deep: bool = False) -> Dict:
        self.flatten_mapping(node)
        keys = set()
        for key_node, _ in node.value:
            key = self.construct_object(key_node, deep=deep)
            try:
                if key in keys:
                    raise yaml.constructor.ConstructorError(
                        "while constructing a mapping", node.start_mark,
                        f"duplicate key: {key!r}", key_node.start_mark,
                    )
                keys.add(key)
            except TypeError as exc:
                raise yaml.constructor.ConstructorError(
                    "while constructing a mapping", node.start_mark,
                    "mapping keys must be hashable", key_node.start_mark,
                ) from exc
        return super().construct_mapping(node, deep=deep)


def _load_yaml_mapping(path: Path) -> Mapping[str, Any]:
    try:
        with path.open(encoding="utf-8") as stream:
            raw = yaml.load(stream, Loader=UniqueKeyLoader)
    except OSError as exc:
        raise ConfigError(f"Cannot read configuration '{path}': {exc}") from exc
    except (yaml.YAMLError, UnicodeError) as exc:
        raise ConfigError(f"Invalid YAML in '{path}': {exc}") from exc

    if not isinstance(raw, Mapping):
        raise ConfigError(f"Configuration '{path}' must contain a YAML mapping")
    return raw


def _format_location(location: Tuple[Union[str, int], ...]) -> str:
    result = ""
    for component in location:
        if isinstance(component, int):
            result += f"[{component}]"
        elif result:
            result += f".{component}"
        else:
            result = component
    return result or "configuration"


def _format_validation_error(path: Path, exc: ValidationError) -> ConfigError:
    messages = []
    for error in exc.errors(
        include_url=False,
        include_context=False,
        include_input=False,
    ):
        location = _format_location(error["loc"])
        message = error["msg"]
        if message.startswith("Value error, "):
            message = message.removeprefix("Value error, ")
        messages.append(f"  - {location}: {message}")
    return ConfigError(
        f"Invalid configuration '{path}':\n" + "\n".join(messages)
    )


def load_case_config(path: Path) -> CaseConfig:
    """Load and validate an experiment against the canonical configuration."""
    raw = _load_yaml_mapping(path)
    try:
        return CaseConfig.model_validate(raw)
    except ValidationError as exc:
        raise _format_validation_error(path, exc) from exc
