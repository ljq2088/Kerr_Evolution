from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import os
from pathlib import Path
import shutil

from kerr_waveform_tools.artifacts import atomic_json
from kerr_waveform_tools.strain import load_ffi_result
from kerr_waveform_tools.waveform_io import load_mode

from .config import EvolutionConfig, OutputRequest
from .errors import ContractError
from .io import load_run_metadata


@dataclass(frozen=True)
class DataRunPaths:
    run_id: str
    staging: Path
    final: Path


def default_data_run_root() -> Path:
    repository = Path(__file__).resolve().parents[4]
    return repository / "data" / "kerr_point_particle_evolution" / "runs"


def make_run_id(config: EvolutionConfig, *, timestamp: datetime | None = None) -> str:
    timestamp = timestamp or datetime.now(timezone.utc)
    spin = f"{config.chi:.12f}".rstrip("0").replace(".", "p")
    return (
        f"kerr_pp_chi{spin}_m{config.m}_nR{config.n_r}_ny{config.n_y}_"
        f"{timestamp.strftime('%Y%m%dT%H%M%SZ')}"
    )


def data_run_paths(root: Path, run_id: str) -> DataRunPaths:
    if not run_id or run_id in (".", "..") or "/" in run_id:
        raise ContractError("run_id must be one nonempty path component")
    root = Path(root)
    return DataRunPaths(run_id, root / ".staging" / run_id, root / run_id)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finalize_data_bundle(
    paths: DataRunPaths,
    config: EvolutionConfig,
    output: OutputRequest,
) -> Path:
    metadata = load_run_metadata(paths.staging)
    if metadata["status"] != "complete" or metadata["run_kind"] != "full":
        raise ContractError("data finalization requires a complete full run")
    products = {}
    for ell in output.ell_out:
        if output.save_psi4_lm:
            mode = paths.staging / f"psi4_l{ell}_m{config.m}.npz"
            if not mode.is_file():
                raise ContractError("requested psi4 mode is missing")
            load_mode(mode, expected_ell=ell, expected_m=config.m)
            figure = paths.staging / f"psi4_l{ell}_m{config.m}_real_abs.png"
            if not figure.is_file() or figure.stat().st_size == 0:
                raise ContractError("requested psi4 plot is missing")
            products[mode.name] = _sha256(mode)
            products[figure.name] = _sha256(figure)
        if output.save_strain_lm:
            strain = paths.staging / f"H_l{ell}_m{config.m}.npz"
            if not strain.is_file():
                raise ContractError("requested strain mode is missing")
            load_ffi_result(strain, expected_ell=ell, expected_m=config.m)
            figure = paths.staging / f"H_l{ell}_m{config.m}_real_abs.png"
            if not figure.is_file() or figure.stat().st_size == 0:
                raise ContractError("requested strain plot is missing")
            products[strain.name] = _sha256(strain)
            products[figure.name] = _sha256(figure)
    if paths.final.exists():
        raise ContractError(f"final data run already exists: {paths.final}")
    checkpoint_directory = paths.staging / "checkpoints"
    if checkpoint_directory.exists():
        shutil.rmtree(checkpoint_directory)
    atomic_json(
        paths.staging / "bundle.json",
        {
            "schema": "kerr-point-particle-data-bundle-v1",
            "status": "complete",
            "review_status": "awaiting_review",
            "run_id": paths.run_id,
            "run_schema": metadata["schema"],
            "config_hash": metadata["config_hash"],
            "code_hash": metadata["code_hash"],
            "trajectory_hash": metadata["trajectory_hash"],
            "products": products,
        },
    )
    paths.final.parent.mkdir(parents=True, exist_ok=True)
    os.replace(paths.staging, paths.final)
    return paths.final


__all__ = [
    "DataRunPaths",
    "data_run_paths",
    "default_data_run_root",
    "finalize_data_bundle",
    "make_run_id",
]
