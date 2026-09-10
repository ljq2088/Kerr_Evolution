from __future__ import annotations

from enum import Enum
import hashlib
import json
import numbers
import os
from pathlib import Path

import numpy as np

from .contracts import ContractError


class ArtifactStatus(str, Enum):
    INCOMPLETE = "incomplete"
    TEST = "test"
    COMPLETE = "complete"


def parse_status(value: str | ArtifactStatus) -> ArtifactStatus:
    try:
        return ArtifactStatus(value)
    except (ValueError, TypeError) as error:
        raise ContractError(f"unknown artifact status: {value!r}") from error


def nonnegative_count(value, name: str) -> int:
    if isinstance(value, np.ndarray):
        if value.shape != ():
            raise ContractError(f"{name} must be a scalar count")
        value = value.item()
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Real):
        raise ContractError(f"{name} must be a non-boolean numeric count")
    numeric = float(value)
    if not np.isfinite(numeric) or numeric < 0 or not numeric.is_integer():
        raise ContractError(f"{name} must be finite, integer-valued, and nonnegative")
    return int(numeric)


def json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(type(value).__name__)


def atomic_json(path: Path, payload: dict) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, sort_keys=True, indent=2, default=json_default) + "\n",
        encoding="ascii",
    )
    os.replace(temporary, path)


def hash_config(payload: str) -> str:
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


_status = parse_status
_nonnegative_count = nonnegative_count
_json_default = json_default
_atomic_json = atomic_json


__all__ = [
    "ArtifactStatus",
    "atomic_json",
    "hash_config",
    "json_default",
    "nonnegative_count",
    "parse_status",
]
