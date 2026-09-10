from __future__ import annotations

from dataclasses import dataclass
import json
import os
from pathlib import Path
import numpy as np

from kerr_waveform_tools.artifacts import (
    ArtifactStatus,
    atomic_json as _atomic_json,
    hash_config,
    nonnegative_count as _nonnegative_count,
    parse_status as _status,
)
from kerr_waveform_tools.waveform_io import MODE_SCHEMA, ModeBuffer, load_mode

from .errors import ContractError
from .evolution import EvolutionState


CHECKPOINT_SCHEMA = "c1-checkpoint-v1"
SCRI_SCHEMA = "c1-scri-v1"
FULL_FIELD_SCHEMA = "c1-full-field-v1"
RUN_SCHEMA = "c1-run-v2"
SUPPORTED_RUN_SCHEMAS = {"c1-run-v1", RUN_SCHEMA}


def _validate_run_counts(metadata: dict) -> dict:
    result = dict(metadata)
    for name in ("N_end", "output_count", "valid_field_slices"):
        if name in result:
            result[name] = _nonnegative_count(result[name], name)
    return result


@dataclass(frozen=True)
class Checkpoint:
    T: float
    state: EvolutionState
    dt: float
    next_output_index: int
    config_hash: str
    trajectory_hash: str
    code_hash: str
    mode_count: int = 0
    scri_count: int = 0
    field_valid_slices: int = 0


def save_checkpoint(path: Path, checkpoint: Checkpoint) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    next_output_index = _nonnegative_count(checkpoint.next_output_index, "next_output_index")
    mode_count = _nonnegative_count(checkpoint.mode_count, "mode_count")
    scri_count = _nonnegative_count(checkpoint.scri_count, "scri_count")
    field_valid_slices = _nonnegative_count(checkpoint.field_valid_slices, "field_valid_slices")
    with temporary.open("wb") as stream:
        np.savez(stream, schema=CHECKPOINT_SCHEMA, T=checkpoint.T, P=checkpoint.state.P, psi=checkpoint.state.psi,
                 status=ArtifactStatus.COMPLETE.value,
                 dt=checkpoint.dt, next_output_index=next_output_index, config_hash=checkpoint.config_hash,
                 trajectory_hash=checkpoint.trajectory_hash, code_hash=checkpoint.code_hash,
                 mode_count=mode_count, scri_count=scri_count, field_valid_slices=field_valid_slices)
    os.replace(temporary, path)


def load_checkpoint(path: Path, *, expected_shape: tuple[int, int], expected_config_hash: str,
                    expected_trajectory_hash: str, expected_code_hash: str) -> Checkpoint:
    with np.load(path, allow_pickle=False) as data:
        if _status(str(data["status"])) is not ArtifactStatus.COMPLETE:
            raise ContractError("checkpoint is not complete")
        if (str(data["schema"]) != CHECKPOINT_SCHEMA or str(data["config_hash"]) != expected_config_hash
                or str(data["trajectory_hash"]) != expected_trajectory_hash or str(data["code_hash"]) != expected_code_hash):
            raise ContractError("checkpoint schema or config hash mismatch")
        P, psi = np.asarray(data["P"]), np.asarray(data["psi"])
        if P.shape != expected_shape or psi.shape != expected_shape or P.dtype != np.complex128 or psi.dtype != np.complex128:
            raise ContractError("checkpoint shape or dtype mismatch")
        return Checkpoint(float(data["T"]), EvolutionState(P, psi), float(data["dt"]),
                          _nonnegative_count(data["next_output_index"], "next_output_index"),
                          expected_config_hash, expected_trajectory_hash, expected_code_hash,
                          _nonnegative_count(data["mode_count"], "mode_count"),
                          _nonnegative_count(data["scri_count"], "scri_count"),
                          _nonnegative_count(data["field_valid_slices"], "field_valid_slices"))


class FullFieldMemmap:
    def __init__(self, path: Path, snapshot_count: int, shape: tuple[int, int]) -> None:
        snapshot_count = _nonnegative_count(snapshot_count, "snapshot_count")
        self.path = Path(path); self.path.parent.mkdir(parents=True, exist_ok=True)
        self.array = np.lib.format.open_memmap(self.path, mode="w+", dtype=np.complex128, shape=(snapshot_count,) + shape)
        self.valid_slices = 0

    def append(self, field: np.ndarray) -> None:
        if self.valid_slices >= len(self.array):
            raise ContractError("full-field memmap capacity exceeded")
        self.array[self.valid_slices] = field; self.valid_slices += 1

    def flush(self) -> None:
        self.array.flush()

    @classmethod
    def resume(cls, path: Path, snapshot_count: int, shape: tuple[int, int], valid_slices: int) -> "FullFieldMemmap":
        snapshot_count = _nonnegative_count(snapshot_count, "snapshot_count")
        valid_slices = _nonnegative_count(valid_slices, "field_valid_slices")
        path = Path(path)
        old = np.load(path, mmap_mode="r")
        if old.shape[1:] != shape or valid_slices > old.shape[0]:
            raise ContractError("full-field restart shape or valid-slice metadata mismatch")
        retained = np.array(old[:valid_slices])
        del old
        temporary = path.with_suffix(path.suffix + ".tmp")
        array = np.lib.format.open_memmap(temporary, mode="w+", dtype=np.complex128,
                                          shape=(max(snapshot_count, valid_slices),) + shape)
        if valid_slices:
            array[:valid_slices] = retained
        array.flush(); del array
        os.replace(temporary, path)
        result = object.__new__(cls)
        result.path = path
        result.array = np.lib.format.open_memmap(path, mode="r+")
        result.valid_slices = valid_slices
        return result


def initialize_run(directory: Path, metadata: dict) -> None:
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=False)
    _atomic_json(directory / "metadata.json", {**_validate_run_counts(metadata), "status": ArtifactStatus.INCOMPLETE.value})


def finalize_run(directory: Path, metadata: dict, *, status: str = "complete") -> None:
    status_value = _status(status)
    if status_value not in (ArtifactStatus.TEST, ArtifactStatus.COMPLETE):
        raise ContractError("final run status must be test or complete")
    _atomic_json(Path(directory) / "metadata.json", {**_validate_run_counts(metadata), "status": status_value.value})


def load_run_metadata(directory: Path, *, allow_incomplete: bool = False) -> dict:
    path = Path(directory) / "metadata.json"
    metadata = json.loads(path.read_text(encoding="ascii"))
    if metadata.get("schema") not in SUPPORTED_RUN_SCHEMAS:
        raise ContractError("run metadata schema mismatch")
    status = _status(metadata.get("status"))
    if status is not ArtifactStatus.COMPLETE and not allow_incomplete:
        raise ContractError("loader refuses an incomplete run")
    return _validate_run_counts(metadata)


def load_scri_field(path: Path, *, expected_config_hash: str | None = None,
                    expected_trajectory_hash: str | None = None, expected_code_hash: str | None = None,
                    allow_incomplete: bool = False) -> tuple[np.ndarray, np.ndarray, dict]:
    with np.load(path, allow_pickle=False) as data:
        T, field = np.asarray(data["T"]), np.asarray(data["psi4_m"])
        metadata = json.loads(str(data["metadata"]))
    if metadata.get("schema") != SCRI_SCHEMA:
        raise ContractError("SCRI field schema mismatch")
    status = _status(metadata.get("status"))
    if status is not ArtifactStatus.COMPLETE and not allow_incomplete:
        raise ContractError("SCRI loader refuses non-complete data")
    for name, expected in (("config_hash", expected_config_hash), ("trajectory_hash", expected_trajectory_hash),
                           ("code_hash", expected_code_hash)):
        if expected is not None and metadata.get(name) != expected:
            raise ContractError(f"SCRI field {name} mismatch")
    if T.ndim != 1 or field.ndim != 2 or field.shape[0] != len(T) or field.dtype != np.complex128:
        raise ContractError("SCRI field shape or dtype mismatch")
    if _nonnegative_count(metadata.get("sample_count"), "sample_count") != len(T):
        raise ContractError("SCRI sample_count does not match stored arrays")
    return T, field, metadata


def write_scri_field(path: Path, T, field, *, metadata: dict, status: str) -> None:
    T, field = np.asarray(T, dtype=float), np.asarray(field, dtype=np.complex128)
    if T.ndim != 1 or field.ndim != 2 or field.shape[0] != len(T):
        raise ContractError("SCRI field write shape mismatch")
    payload = {**metadata, "schema": SCRI_SCHEMA, "status": _status(status).value, "sample_count": len(T)}
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, T=T, psi4_m=field, metadata=json.dumps(payload, sort_keys=True))
    os.replace(temporary, path)


def write_full_field_metadata(path: Path, *, status: str, valid_slices: int, shape: tuple[int, ...],
                              config_hash: str, trajectory_hash: str, code_hash: str) -> None:
    valid_slices = _nonnegative_count(valid_slices, "valid_slices")
    checked_shape = [_nonnegative_count(value, "shape entry") for value in shape]
    _atomic_json(Path(path), {"schema": FULL_FIELD_SCHEMA, "status": _status(status).value, "valid_slices": valid_slices,
                             "shape": checked_shape, "dtype": "complex128", "config_hash": config_hash,
                             "trajectory_hash": trajectory_hash, "code_hash": code_hash})


def load_full_field(path: Path, metadata_path: Path, *, expected_config_hash: str | None = None,
                    expected_trajectory_hash: str | None = None, expected_code_hash: str | None = None,
                    allow_incomplete: bool = False) -> tuple[np.ndarray, dict]:
    metadata = json.loads(Path(metadata_path).read_text(encoding="ascii"))
    if metadata.get("schema") != FULL_FIELD_SCHEMA or metadata.get("dtype") != "complex128":
        raise ContractError("full-field metadata schema or dtype mismatch")
    status = _status(metadata.get("status"))
    if status is not ArtifactStatus.COMPLETE and not allow_incomplete:
        raise ContractError("full-field loader refuses non-complete data")
    for name, expected in (("config_hash", expected_config_hash), ("trajectory_hash", expected_trajectory_hash),
                           ("code_hash", expected_code_hash)):
        if expected is not None and metadata.get(name) != expected:
            raise ContractError(f"full-field {name} mismatch")
    array = np.load(path, mmap_mode="r")
    valid = _nonnegative_count(metadata.get("valid_slices"), "valid_slices")
    stored_shape = tuple(_nonnegative_count(value, "shape entry") for value in metadata.get("shape", ()))
    if stored_shape != array.shape or array.dtype != np.complex128 or valid > array.shape[0]:
        raise ContractError("full-field shape or valid-slice metadata mismatch")
    return array[:valid], metadata


class RollingCheckpoint:
    def __init__(self, directory: Path) -> None:
        self.directory = Path(directory); self.directory.mkdir(parents=True, exist_ok=True)

    def save(self, checkpoint: Checkpoint) -> Path:
        marker = self.directory / "latest.json"
        generation = 0
        if marker.exists():
            generation = _nonnegative_count(json.loads(marker.read_text(encoding="ascii"))["generation"], "generation") + 1
        path = self.directory / f"checkpoint_{generation % 2}.npz"
        save_checkpoint(path, checkpoint)
        _atomic_json(marker, {"generation": generation, "path": path.name})
        return path

    def latest(self, **expected) -> Checkpoint:
        marker = self.directory / "latest.json"
        if not marker.exists():
            raise ContractError("no rolling checkpoint is available")
        name = json.loads(marker.read_text(encoding="ascii"))["path"]
        return load_checkpoint(self.directory / name, **expected)
