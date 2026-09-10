from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

from .artifacts import ArtifactStatus, json_default, nonnegative_count, parse_status
from .contracts import ContractError


MODE_SCHEMA = "c1-mode-v2"
SUPPORTED_MODE_SCHEMAS = {"c1-mode-v1", MODE_SCHEMA}


class ModeBuffer:
    REQUIRED_METADATA = {
        "ell",
        "m",
        "events",
        "normalization",
        "error_metadata",
        "input_hash",
        "code_hash",
        "config_hash",
        "analysis_start",
        "phi_dot_at_analysis_start",
        "run_kind",
        "T_bound",
        "T1",
        "N_end",
        "delta_t_post",
        "environment",
    }

    def __init__(self, path: Path, metadata: dict, chunk_size: int = 1024) -> None:
        chunk_size = nonnegative_count(chunk_size, "chunk_size")
        if chunk_size == 0:
            raise ContractError("mode chunk size must be positive")
        missing = self.REQUIRED_METADATA - metadata.keys()
        if missing:
            raise ContractError(f"mode metadata missing keys: {sorted(missing)}")
        metadata = dict(metadata)
        metadata["N_end"] = nonnegative_count(metadata["N_end"], "N_end")
        self.path, self.chunk_size = Path(path), chunk_size
        self.metadata = dict(metadata)
        self.times: list[float] = []
        self.values: list[complex] = []

    def append(self, T: float, value: complex) -> None:
        self.times.append(float(T))
        self.values.append(complex(value))
        if len(self.times) % self.chunk_size == 0:
            self.flush(complete=False)

    def flush(self, *, complete: bool = False, status: str | None = None) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.path.with_suffix(self.path.suffix + ".tmp")
        status_value = parse_status(
            status or (ArtifactStatus.COMPLETE.value if complete else ArtifactStatus.INCOMPLETE.value)
        ).value
        metadata = {**self.metadata, "schema": MODE_SCHEMA, "status": status_value}
        metadata["sample_count"] = len(self.times)
        with temporary.open("wb") as stream:
            np.savez_compressed(
                stream,
                T=np.asarray(self.times),
                psi4_lm=np.asarray(self.values, dtype=np.complex128),
                metadata=json.dumps(metadata, sort_keys=True, default=json_default),
            )
        os.replace(temporary, self.path)

    @classmethod
    def resume(
        cls,
        path: Path,
        metadata: dict,
        chunk_size: int,
        *,
        truncate: int | None = None,
    ) -> "ModeBuffer":
        T, values, stored = load_mode(path, allow_incomplete=True)
        for key in cls.REQUIRED_METADATA:
            if stored[key] != metadata[key]:
                raise ContractError(f"mode restart metadata mismatch for {key}")
        result = cls(path, metadata, chunk_size)
        if truncate is not None:
            truncate = nonnegative_count(truncate, "mode_count")
            if len(T) < truncate:
                raise ContractError("mode file contains fewer samples than checkpoint mode_count")
            T, values = T[:truncate], values[:truncate]
        result.times = T.tolist()
        result.values = values.tolist()
        if truncate is not None:
            result.flush(status=ArtifactStatus.INCOMPLETE.value)
            verified_T, _, _ = load_mode(path, allow_incomplete=True)
            if len(verified_T) != truncate:
                raise ContractError("atomic mode truncation did not reach checkpoint mode_count")
        return result


def load_mode(
    path: Path,
    *,
    expected_ell: int | None = None,
    expected_m: int | None = None,
    expected_environment: dict | None = None,
    allow_incomplete: bool = False,
) -> tuple[np.ndarray, np.ndarray, dict]:
    with np.load(path, allow_pickle=False) as data:
        T, values = np.asarray(data["T"]), np.asarray(data["psi4_lm"])
        metadata = json.loads(str(data["metadata"]))
    missing = ModeBuffer.REQUIRED_METADATA - metadata.keys()
    if metadata.get("schema") not in SUPPORTED_MODE_SCHEMAS or missing:
        raise ContractError(f"mode schema invalid or metadata missing: {sorted(missing)}")
    status = parse_status(metadata.get("status"))
    if status is not ArtifactStatus.COMPLETE and not allow_incomplete:
        raise ContractError("mode loader refuses an incomplete run")
    if (
        (expected_ell is not None and metadata["ell"] != expected_ell)
        or (expected_m is not None and metadata["m"] != expected_m)
    ):
        raise ContractError("mode indices do not match request")
    if expected_environment is not None and metadata["environment"] != expected_environment:
        raise ContractError("mode environment mismatch")
    if T.ndim != 1 or values.shape != T.shape or values.dtype != np.complex128:
        raise ContractError("mode file shape or dtype mismatch")
    if nonnegative_count(metadata.get("sample_count"), "sample_count") != len(T):
        raise ContractError("mode sample_count does not match stored arrays")
    metadata["N_end"] = nonnegative_count(metadata["N_end"], "N_end")
    return T, values, metadata


__all__ = ["MODE_SCHEMA", "SUPPORTED_MODE_SCHEMAS", "ModeBuffer", "load_mode"]
