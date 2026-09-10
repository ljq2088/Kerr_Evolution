from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import numpy as np
from dataclasses import asdict
import sys
import scipy
from tqdm.auto import tqdm
import kerr_waveform_tools
from kerr_waveform_tools.plotting import save_complex_series_plot
from kerr_waveform_tools.strain import (
    fixed_frequency_integrate,
    load_ffi_result,
    write_ffi_result,
)
from kerr_waveform_tools.waveform_io import load_mode

from .config import EvolutionConfig, OutputRequest
from .evolution import EvolutionState, SpatialRHS, TimeStepChoice, choose_timestep, rk4_step, zero_state
from .extraction import project_scri_mode, scri_slice
from .grid import SpatialGrid, build_spatial_grid
from .initial_data import PlungeInitialData, validate_initial_data
from .io import (Checkpoint, FullFieldMemmap, ModeBuffer, RollingCheckpoint, finalize_run,
                 RUN_SCHEMA, hash_config, initialize_run, load_full_field, load_scri_field, write_full_field_metadata,
                 write_scri_field)
from .source import GridSourceEvaluator
from .trajectory import DOP853Trajectory, integrate_trajectory


@dataclass(frozen=True)
class SingleSystem:
    config: EvolutionConfig
    output: OutputRequest
    initial_data: PlungeInitialData
    grid: SpatialGrid
    trajectory: DOP853Trajectory
    source: GridSourceEvaluator
    rhs: SpatialRHS
    timestep: TimeStepChoice
    config_hash: str
    trajectory_hash: str
    code_hash: str
    T_bound: float
    T1: float
    N_end: int
    initial_data_provenance: dict

    @property
    def initial_state(self):
        return zero_state(self.grid.RR.shape)


def assemble_single_system(
    config: EvolutionConfig,
    initial_data: PlungeInitialData,
    output: OutputRequest | None = None,
    *,
    initial_data_provenance: dict | None = None,
) -> SingleSystem:
    output = output or OutputRequest(ell_out=(config.m,))
    if any(ell < abs(config.m) for ell in output.ell_out):
        raise ValueError("each requested ell must satisfy ell >= |m|")
    validate_initial_data(config, initial_data)
    grid = build_spatial_grid(config.chi, config.m, config.n_r, config.n_y)
    sigma_R, sigma_y = config.sigma_r_over_dr * grid.dR, config.sigma_y_over_dy * grid.dy
    trajectory = integrate_trajectory(initial_data, config, sigma_R)
    source = GridSourceEvaluator(grid, trajectory, config.m, sigma_R, sigma_y,
                                 config.gaussian_tail_tolerance, config.tau_on_over_m)
    rhs = SpatialRHS(grid, config.chi, config.m, source)
    timestep = choose_timestep(grid, config.chi, config.m, trajectory, output.output_dt_over_m,
                               tau_on=config.tau_on_over_m, sigma_R=sigma_R)
    config_hash = hash_config(config.to_json() + output.to_json())
    trajectory_hash = hashlib.sha256(np.asarray(trajectory.T).tobytes() + np.asarray(trajectory.R).tobytes()
                                     + np.asarray(trajectory.Phi).tobytes()).hexdigest()
    T_bound = trajectory.evolution_end_time(config.post_light_ring_time_over_m)
    N_end = int(np.floor(T_bound / output.output_dt_over_m))
    T1 = N_end * output.output_dt_over_m
    return SingleSystem(config, output, initial_data, grid, trajectory, source, rhs, timestep,
                        config_hash, trajectory_hash, _source_code_hash(), T_bound, T1, N_end,
                        dict(initial_data_provenance or {}))


def _source_code_hash() -> str:
    digest = hashlib.sha256()
    packages = (
        Path(__file__).resolve().parent,
        Path(kerr_waveform_tools.__file__).resolve().parent,
    )
    for package in packages:
        for path in sorted(package.glob("*.py")):
            digest.update(package.name.encode("ascii"))
            digest.update(path.name.encode("ascii"))
            digest.update(path.read_bytes())
    return digest.hexdigest()


@dataclass(frozen=True)
class RunResult:
    directory: Path
    final_time: float
    final_state: EvolutionState
    output_count: int


def _mode_metadata(system: SingleSystem, ell: int, run_kind: str) -> dict:
    input_hash = hashlib.sha256(json.dumps(system.initial_data.as_record(), sort_keys=True).encode("ascii")).hexdigest()
    return {
        "ell": ell, "m": system.config.m, "events": asdict(system.trajectory.events),
        "normalization": "Ripley peeling psi4 at SCRI+, 2pi spin-weighted projection",
        "error_metadata": {"V1": "open", "V2": "not_run"}, "input_hash": input_hash,
        "code_hash": system.code_hash, "config_hash": system.config_hash,
        "run_kind": run_kind, "T_bound": system.T_bound, "T1": system.T1,
        "N_end": system.N_end, "delta_t_post": system.config.post_light_ring_time_over_m,
        "environment": _environment(),
        "initial_data_provenance": system.initial_data_provenance,
        "analysis_start": system.timestep.t_stable,
        "phi_dot_at_analysis_start": system.trajectory.sample(system.timestep.t_stable).Phi_jet.first,
        "light_ring_scri_reference": system.trajectory.light_ring_scri_reference,
        "background": {
            "mass_scale": system.config.mass_scale,
            "M_internal": system.config.M,
            "chi": system.config.chi,
            "a_internal": system.config.a,
            "L_internal": system.config.L,
        },
        "amplitude_scaling": (
            "unit particle-mass response; no physical mu, luminosity-distance, "
            "or BBH amplitude mapping applied"
        ),
    }


def _environment() -> dict:
    return {"python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__}


def finalize_mode_products(system: SingleSystem, directory: Path, *, full_run: bool) -> None:
    if not full_run:
        return
    directory = Path(directory)
    for ell in system.output.ell_out:
        mode_path = directory / f"psi4_l{ell}_m{system.config.m}.npz"
        if not mode_path.exists():
            continue
        mode_T, mode_values, mode_metadata = load_mode(
            mode_path,
            expected_ell=ell,
            expected_m=system.config.m,
        )
        if system.output.save_psi4_lm:
            save_complex_series_plot(
                mode_T,
                mode_values,
                directory / f"psi4_l{ell}_m{system.config.m}_real_abs.png",
                symbol=rf"\psi_{{4,{ell}{system.config.m}}}",
                title=rf"Peeling waveform: $\chi={system.config.chi:g}$, "
                      rf"$(\ell,m)=({ell},{system.config.m})$",
                event_lines={
                    r"$T_{\rm stable}$": system.timestep.t_stable,
                    r"$T_{\rm LR}^{\mathcal{I}^+}$": system.trajectory.light_ring_scri_reference,
                    r"$T_{\rm source\ off}$": system.trajectory.events.source_off,
                },
            )
        if system.output.save_strain_lm:
            ffi = fixed_frequency_integrate(
                mode_T,
                mode_values,
                ell=ell,
                m=system.config.m,
                phi_dot_at_start=mode_metadata["phi_dot_at_analysis_start"],
                transform_start=float(mode_T[0]),
                analysis_start=mode_metadata["analysis_start"],
                stop=mode_metadata["T1"],
                input_schema=mode_metadata["schema"],
                source_input_hash=mode_metadata["input_hash"],
                source_events=mode_metadata["events"],
                source_normalization=mode_metadata["normalization"],
                source_error_metadata=mode_metadata["error_metadata"],
                source_environment=mode_metadata["environment"],
                source_background=mode_metadata["background"],
                amplitude_scaling=mode_metadata["amplitude_scaling"],
                sampling_status="owner-waived",
            )
            strain_path = directory / f"H_l{ell}_m{system.config.m}.npz"
            write_ffi_result(
                strain_path,
                ffi,
                status="complete",
                source_mode_metadata=mode_metadata,
            )
            loaded_ffi, strain_metadata = load_ffi_result(
                strain_path,
                expected_ell=ell,
                expected_m=system.config.m,
            )
            save_complex_series_plot(
                loaded_ffi.T,
                loaded_ffi.H,
                directory / f"H_l{ell}_m{system.config.m}_real_abs.png",
                symbol=rf"H_{{{ell}{system.config.m}}}",
                title=rf"FFI strain: $\chi={system.config.chi:g}$, "
                      rf"$(\ell,m)=({ell},{system.config.m})$",
                event_lines={
                    r"$T_{\rm stable}$": strain_metadata["analysis_start"],
                    r"$T_{\rm LR}^{\mathcal{I}^+}$": system.trajectory.light_ring_scri_reference,
                },
                shaded_intervals={
                    "FFI taper": (
                        strain_metadata["transform_start"],
                        strain_metadata["transform_start"] + strain_metadata["taper_width"],
                    ),
                    "FFI end taper": (
                        strain_metadata["recommended_stop"],
                        strain_metadata["transform_stop"],
                    ),
                },
            )
            if not system.output.save_psi4_lm:
                mode_path.unlink()


def run_single_system(system: SingleSystem, directory: Path, *, T_end: float | None = None,
                      restart: bool = False, show_progress: bool = True) -> RunResult:
    directory = Path(directory)
    T_end = system.T1 if T_end is None else float(T_end)
    if T_end <= 0:
        raise ValueError("T_end must be positive")
    if T_end > system.T_bound + 1e-13:
        raise ValueError("T_end cannot exceed the approved post-light-ring upper bound")
    full_run = np.isclose(T_end, system.T1, rtol=0.0, atol=1e-13)
    if system.output.save_strain_lm and not full_run:
        raise ValueError("strain output requires a complete run through T1")
    metadata = {
        "schema": RUN_SCHEMA, "config": json.loads(system.config.to_json()),
        "output": json.loads(system.output.to_json()), "initial_data": system.initial_data.as_record(),
        "events": asdict(system.trajectory.events), "light_ring_scri_reference": system.trajectory.light_ring_scri_reference,
        "T_bound": system.T_bound, "T1": system.T1, "N_end": system.N_end,
        "delta_t_post": system.config.post_light_ring_time_over_m,
        "T_end": T_end, "run_kind": "full" if full_run else "test",
        "grid_version": system.grid.version, "stencil_width": 7,
        "timestep": asdict(system.timestep), "config_hash": system.config_hash,
        "trajectory_hash": system.trajectory_hash, "code_hash": system.code_hash,
        "environment": _environment(),
        "initial_data_provenance": system.initial_data_provenance,
    }
    shape = system.grid.RR.shape
    if restart:
        rolling = RollingCheckpoint(directory / "checkpoints")
        existing = json.loads((directory / "metadata.json").read_text(encoding="ascii"))
        for key in ("schema", "config_hash", "trajectory_hash", "code_hash", "environment"):
            if existing.get(key) != metadata.get(key):
                raise ValueError(f"restart metadata mismatch for {key}")
        checkpoint = rolling.latest(expected_shape=shape, expected_config_hash=system.config_hash,
                                    expected_trajectory_hash=system.trajectory_hash, expected_code_hash=system.code_hash)
        state, T, output_index = checkpoint.state, checkpoint.T, checkpoint.next_output_index
    else:
        initialize_run(directory, metadata)
        state, T, output_index = system.initial_state, 0.0, 0
        rolling = RollingCheckpoint(directory / "checkpoints")
    mode_buffers = {}
    if system.output.save_psi4_lm or system.output.save_strain_lm:
        for ell in system.output.ell_out:
            path = directory / f"psi4_l{ell}_m{system.config.m}.npz"
            mode_meta = _mode_metadata(system, ell, "full" if full_run else "test")
            if restart and output_index and not path.exists():
                raise ValueError("restart checkpoint has outputs but mode file is missing")
            mode_buffers[ell] = (ModeBuffer.resume(path, mode_meta, system.output.chunk_size, truncate=checkpoint.mode_count)
                                 if restart and path.exists() else ModeBuffer(path, mode_meta, system.output.chunk_size))
    scri_times, scri_values = [], []
    if restart and system.output.save_scri_field and (directory / "scri_field.npz").exists():
        stored_T, stored_field, _ = load_scri_field(
            directory / "scri_field.npz", expected_config_hash=system.config_hash,
            expected_trajectory_hash=system.trajectory_hash, expected_code_hash=system.code_hash,
            allow_incomplete=True)
        scri_times = stored_T[:checkpoint.scri_count].tolist()
        scri_values = stored_field[:checkpoint.scri_count].tolist()
        if len(stored_T) < checkpoint.scri_count:
            raise ValueError("SCRI file contains fewer samples than checkpoint scri_count")
        if len(stored_T) > checkpoint.scri_count:
            write_scri_field(directory / "scri_field.npz", scri_times, scri_values,
                             metadata={"config_hash": system.config_hash,
                                       "trajectory_hash": system.trajectory_hash,
                                       "code_hash": system.code_hash}, status="incomplete")
            verified_T, _, _ = load_scri_field(directory / "scri_field.npz", allow_incomplete=True)
            if len(verified_T) != checkpoint.scri_count:
                raise ValueError("atomic SCRI truncation did not reach checkpoint scri_count")
    elif restart and system.output.save_scri_field and checkpoint.scri_count:
        raise ValueError("restart checkpoint has SCRI samples but SCRI file is missing")
    field_writer = None
    if system.output.save_field:
        count = max(1, int(np.floor(T_end / system.output.field_dt_over_m + 1e-12)))
        if restart:
            available_field, _ = load_full_field(
                directory / "field.npy", directory / "field_metadata.json",
                expected_config_hash=system.config_hash, expected_trajectory_hash=system.trajectory_hash,
                expected_code_hash=system.code_hash, allow_incomplete=True)
            if len(available_field) < checkpoint.field_valid_slices:
                raise ValueError("full-field file contains fewer slices than checkpoint field_valid_slices")
        field_writer = (FullFieldMemmap.resume(directory / "field.npy", count, shape,
                                               checkpoint.field_valid_slices)
                        if restart else FullFieldMemmap(directory / "field.npy", count, shape))
        if restart:
            write_full_field_metadata(directory / "field_metadata.json", status="incomplete",
                                      valid_slices=field_writer.valid_slices, shape=field_writer.array.shape,
                                      config_hash=system.config_hash, trajectory_hash=system.trajectory_hash,
                                      code_hash=system.code_hash)
            verified_field, _ = load_full_field(directory / "field.npy", directory / "field_metadata.json",
                                                allow_incomplete=True)
            if len(verified_field) != checkpoint.field_valid_slices:
                raise ValueError("atomic full-field truncation did not reach checkpoint field_valid_slices")
    field_every = round(system.output.field_dt_over_m / system.output.output_dt_over_m)
    checkpoint_every = round(system.output.checkpoint_dt_over_m / system.output.output_dt_over_m)
    next_output_time = (output_index + 1) * system.output.output_dt_over_m

    def flush_scri(status: str) -> None:
        if not system.output.save_scri_field:
            return
        write_scri_field(directory / "scri_field.npz", scri_times, scri_values,
                         metadata={"config_hash": system.config_hash, "trajectory_hash": system.trajectory_hash,
                                   "code_hash": system.code_hash}, status=status)

    def save_restart_point() -> None:
        for buffer in mode_buffers.values():
            buffer.flush(complete=False)
        flush_scri("incomplete")
        if field_writer is not None:
            field_writer.flush()
            write_full_field_metadata(directory / "field_metadata.json", status="incomplete",
                                      valid_slices=field_writer.valid_slices, shape=field_writer.array.shape,
                                      config_hash=system.config_hash, trajectory_hash=system.trajectory_hash,
                                      code_hash=system.code_hash)
        rolling.save(Checkpoint(T, state, system.timestep.dt, output_index, system.config_hash,
                                system.trajectory_hash, system.code_hash,
                                min((len(buffer.times) for buffer in mode_buffers.values()), default=0),
                                len(scri_times), 0 if field_writer is None else field_writer.valid_slices))

    rk_steps = 0
    with tqdm(total=T_end, initial=T, disable=not show_progress, desc=f"m={system.config.m} RK4",
              unit="M", dynamic_ncols=True, mininterval=0.5) as progress:
        while T < T_end and not np.isclose(T, T_end, rtol=0.0, atol=1e-15):
            until_output = next_output_time - T
            dt = min(system.timestep.dt, T_end - T, until_output if until_output > 1e-15 else system.timestep.dt)
            state = rk4_step(system.rhs, T, state, dt)
            T += dt
            rk_steps += 1
            progress.update(dt)
            if np.isclose(T, next_output_time, rtol=0.0, atol=2e-13):
                T = next_output_time
                progress.set_postfix(T=f"{T:.2f}", steps=rk_steps, refresh=False)
                if system.output.save_scri_field:
                    scri_times.append(T); scri_values.append(scri_slice(state.psi, system.grid.R))
                for ell, buffer in mode_buffers.items():
                    buffer.append(T, project_scri_mode(state.psi, system.grid, ell, system.config.m))
                output_index += 1
                next_output_time = (output_index + 1) * system.output.output_dt_over_m
                if field_writer is not None and output_index % field_every == 0:
                    field_writer.append(state.psi); field_writer.flush()
                if output_index % checkpoint_every == 0:
                    save_restart_point()
    save_restart_point()
    for buffer in mode_buffers.values():
        buffer.flush(status="complete" if full_run else "test")
    finalize_mode_products(system, directory, full_run=full_run)
    if system.output.save_scri_field:
        flush_scri("complete" if full_run else "test")
    if field_writer is not None:
        field_writer.flush()
        metadata["valid_field_slices"] = field_writer.valid_slices
        write_full_field_metadata(directory / "field_metadata.json",
                                  status="complete" if full_run else "test",
                                  valid_slices=field_writer.valid_slices, shape=field_writer.array.shape,
                                  config_hash=system.config_hash, trajectory_hash=system.trajectory_hash,
                                  code_hash=system.code_hash)
    metadata["output_count"] = output_index
    metadata["final_time"] = T
    metadata["final_time_on_output_grid"] = bool(np.isclose(T / system.output.output_dt_over_m,
                                                            round(T / system.output.output_dt_over_m), atol=1e-12))
    finalize_run(directory, metadata, status="complete" if full_run else "test")
    return RunResult(directory, T, state, output_index)
