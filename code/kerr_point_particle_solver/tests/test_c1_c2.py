import json
import math
from pathlib import Path
from dataclasses import asdict, replace
import subprocess
import sys
import hashlib
import shutil

import numpy as np
import pytest
from scipy import sparse
import scipy
import kerr_waveform_tools

from sminus2_point_particle.config import EvolutionConfig, OutputRequest
from sminus2_point_particle.evolution import EvolutionState, SpatialRHS, choose_timestep, evolve_accepted_steps, rk4_step, zero_state
from sminus2_point_particle.extraction import project_scri_mode
from sminus2_point_particle.grid import build_spatial_grid, reconstruct_axis_endpoints, simpson_weights_augmented
from sminus2_point_particle.initial_data import FluxRecord, JsonInitialDataCache, PlungeInitialData, generate_ori_thorne_initial_data, validate_initial_data
from sminus2_point_particle.io import (
    Checkpoint,
    FullFieldMemmap,
    ModeBuffer,
    RollingCheckpoint,
    finalize_run,
    hash_config,
    load_checkpoint,
    load_full_field,
    load_mode,
    load_run_metadata,
    load_scri_field,
    save_checkpoint,
)
from sminus2_point_particle.operators import finite_difference_weights
from sminus2_point_particle.source import (
    GridSourceEvaluator,
    evaluate_regrouped_horizon_total_limit,
    evaluate_source_blocks,
    quintic_window,
    window_amplitude,
)
from sminus2_point_particle.run_storage import data_run_paths, finalize_data_bundle
from sminus2_point_particle.strain import absolute_coarse_indices, compare_sampling, fixed_frequency_integrate, load_ffi_result, quintic_taper
from sminus2_point_particle.trajectory import Jet3, integrate_trajectory
from sminus2_point_particle.workflow import (
    assemble_single_system,
    finalize_mode_products,
    run_single_system,
)


FIXTURE = Path(__file__).parent / "fixtures" / "chi0p8_l6_flux.json"


@pytest.fixture(scope="module")
def plunge():
    return generate_ori_thorne_initial_data(FluxRecord(**json.loads(FIXTURE.read_text(encoding="ascii"))))


@pytest.fixture(scope="module")
def grid():
    return build_spatial_grid(0.8, 2, 128, 33)


@pytest.fixture(scope="module")
def trajectory(plunge, grid):
    config = EvolutionConfig(1, 0.8, 2, 128, 33, flux_model_version=plunge.manifest.flux_model_version)
    return integrate_trajectory(plunge, config, 4 * grid.dR)


def test_target_consistent_initial_data_fixture_and_atomic_cache(tmp_path, plunge):
    assert plunge.r0 == pytest.approx(2.894197755074886)
    assert plunge.energy == pytest.approx(0.877792349241831)
    assert plunge.angular_momentum == pytest.approx(2.3800439756397718)
    manifest = plunge.manifest
    assert (manifest.q_mass, manifest.lmax, manifest.t_ot, manifest.x_start) == (1e-5, 6, 3.412, -0.5381356808642278)
    cfg = EvolutionConfig(1, 0.8, 2, 128, 33, flux_model_version=manifest.flux_model_version)
    key = cfg.initial_data_key(manifest.schema_version, manifest.generator_version)
    cache = JsonInitialDataCache(tmp_path, cfg)
    assert cache.load(key) is None
    cache.store(key, plunge)
    assert cache.load(key) == plunge
    cache.store(key, plunge)
    with pytest.raises(Exception):
        validate_initial_data(cfg, PlungeInitialData(plunge.r0, plunge.energy, plunge.angular_momentum,
                                                    replace(manifest, alpha_ot=manifest.alpha_ot * 1.01)))
    with pytest.raises(Exception):
        cache.store("0" * 64, plunge)
    record_path = tmp_path / f"{key}.json"
    record = json.loads(record_path.read_text(encoding="ascii")); record["manifest"]["chi"] = 0.7
    record_path.write_text(json.dumps(record), encoding="ascii")
    with pytest.raises(Exception):
        cache.load(key)


def test_dop853_unique_trajectory_events_and_regular_horizon_extension(plunge, grid, trajectory):
    events = trajectory.events
    assert events.light_ring < events.horizon_crossing == events.source_off < events.terminal
    assert trajectory.evolution_end_time() == pytest.approx(trajectory.light_ring_scri_reference + 120)
    before = trajectory.sample(events.horizon_crossing - 1e-6)
    after = trajectory.sample(events.horizon_crossing + 1e-6)
    for left, right in ((before.R_jet, after.R_jet), (before.Phi_jet, after.Phi_jet), (before.u_T, after.u_T)):
        assert np.all(np.isfinite([left.value, left.first, left.second, right.value, right.first, right.second]))
        assert right.value == pytest.approx(left.value, rel=2e-6, abs=2e-6)
    config = EvolutionConfig(1, 0.8, 2, 128, 33, flux_model_version=plunge.manifest.flux_model_version)
    refined = integrate_trajectory(plunge, config, 4 * grid.dR, rtol=2e-12, atol=2e-14, max_step=0.125)
    assert refined.events.light_ring == pytest.approx(events.light_ring, abs=2e-8)
    assert refined.events.horizon_crossing == pytest.approx(events.horizon_crossing, abs=2e-8)
    with pytest.raises(Exception):
        integrate_trajectory(plunge, replace(config, chi=0.7), 4 * grid.dR)


def test_seven_point_grid_polynomials_factor_rows_and_simpson(grid):
    for degree in range(7):
        f = grid.R**degree
        d1 = np.zeros_like(f) if degree == 0 else degree * grid.R ** (degree - 1)
        d2 = np.zeros_like(f) if degree < 2 else degree * (degree - 1) * grid.R ** (degree - 2)
        assert np.allclose(grid.D1R @ f, d1, atol=2e-9)
        assert np.allclose(grid.D2R @ f, d2, atol=2e-6)
    f = (1 - grid.y) ** grid.p_south * (1 + 0.3 * grid.y + 0.2 * grid.y**2)
    f1 = -(grid.p_south) * (1 - grid.y) ** (grid.p_south - 1) * (1 + 0.3 * grid.y + 0.2 * grid.y**2) + (1 - grid.y) ** grid.p_south * (0.3 + 0.4 * grid.y)
    assert np.allclose((grid.D1y_field @ f)[-3:], f1[-3:], atol=2e-11)
    augmented = reconstruct_axis_endpoints(f, grid)
    assert augmented[-1] == 0
    nodes, weights = simpson_weights_augmented(grid)
    assert np.sum(weights) == pytest.approx(2.0)
    assert abs(np.sum(weights * nodes**4) - 2 / 5) < 4e-6


def test_sparse_scaled_operators_m2_m4_and_medium_construction():
    conditions = {}
    for m in (2, 4):
        small = build_spatial_grid(0.8, m, 128, 33)
        medium = build_spatial_grid(0.8, m, 512, 129)
        for candidate in (small, medium):
            assert all(sparse.isspmatrix_csr(matrix) for matrix in
                       (candidate.D1R, candidate.D2R, candidate.D1y_field, candidate.D2y_field,
                        candidate.D1y_source, candidate.D2y_source))
            assert candidate.D1R.nnz <= 7 * candidate.R.size
            assert candidate.D1y_field.nnz <= 7 * candidate.y.size
            assert candidate.max_factor_condition < 1e6
        conditions[m] = (small.max_factor_condition, medium.max_factor_condition)
        assert conditions[m][1] == pytest.approx(conditions[m][0], rel=1e-9)
        p_n, p_s = small.p_north, small.p_south
        field = (1 + small.y) ** p_n * (1 - small.y) ** p_s
        assert np.all(np.isfinite(small.D1y_field @ field))


def test_startup_cutoff_horizon_limit_and_source_off(grid, trajectory):
    window = quintic_window(10.0)
    assert window.value == pytest.approx(0.5)
    assert quintic_window(0) == Jet3(0, 0, 0)
    assert quintic_window(20) == Jet3(1, 0, 0)
    amplitude = Jet3(2 + 1j, 3 - 2j, -1 + 4j)
    dressed = window_amplitude(amplitude, window)
    assert dressed.first == pytest.approx(window.first * amplitude.value + window.value * amplitude.first)
    source = GridSourceEvaluator(grid, trajectory, 2, 4 * grid.dR, 4 * grid.dy)
    assert np.max(np.abs(source.blocks(0).total)) == 0
    exterior = source.blocks(trajectory.events.horizon_crossing - 1e-4).total
    assert np.all(np.isfinite(exterior)) and np.max(np.abs(exterior[-1])) == 0
    cols = np.arange(1, 8)
    extrapolation = np.tensordot(finite_difference_weights(grid.R[0], grid.R[cols], 0), exterior[cols], axes=(0, 0))
    assert np.linalg.norm(exterior[0] - extrapolation) / np.linalg.norm(exterior[0]) < 0.05
    assert np.array_equal(source.blocks(trajectory.events.source_off).total, np.zeros_like(exterior))


@pytest.mark.parametrize("chi,m,flux_value", [(0.2, 2, 0.0015), (0.95, 4, 0.08)])
def test_regrouped_horizon_limit_multiple_spins_modes_resolutions(chi, m, flux_value):
    initial = generate_ori_thorne_initial_data(FluxRecord(chi, 6, flux_value, 0.01, f"synthetic-{chi}",
                                                          "synthetic numerical-contract fixture only"))
    for n_r in (128, 192):
        candidate = build_spatial_grid(chi, m, n_r, 33)
        config = EvolutionConfig(1, chi, m, n_r, 33, flux_model_version=initial.manifest.flux_model_version)
        plunge_trajectory = integrate_trajectory(initial, config, 4 * candidate.dR)
        source = GridSourceEvaluator(candidate, plunge_trajectory, m, 4 * candidate.dR, 4 * candidate.dy)
        for offset in (0.1, 0.01, 0.001):
            T = plunge_trajectory.events.horizon_crossing - offset
            blocks = source.blocks(T)
            reference = evaluate_source_blocks(
                candidate.RR, candidate.yy, plunge_trajectory.sample(T), 1.0, chi, 1.0, m,
                source.sigma_R, source.sigma_y, candidate.dr, candidate.dy_source,
                epsilon_g=source.epsilon_g, allow_boundary_limits=True,
            )
            for actual, expected in zip(
                (blocks.block1, blocks.block2, blocks.block3, blocks.block4, blocks.total),
                (reference.block1, reference.block2, reference.block3, reference.block4, reference.total),
            ):
                assert np.allclose(actual, expected, rtol=2e-12, atol=2e-12)
            assert blocks.horizon_total_limit is not None
            assert np.all(np.isfinite(blocks.total[0]))
            assert np.array_equal(blocks.total[0], blocks.horizon_total_limit[0])
            independent = evaluate_regrouped_horizon_total_limit(
                candidate.RR, candidate.yy, plunge_trajectory.sample(T), 1.0, chi, 1.0, m,
                source.sigma_R, source.sigma_y, candidate.dr, candidate.dy_source,
                source.epsilon_g,
            )
            assert np.allclose(blocks.horizon_total_limit[0], independent[0], rtol=2e-12, atol=2e-12)
            generic_sum = blocks.block1[0] + blocks.block2[0] + blocks.block3[0] + blocks.block4[0]
            assert np.allclose(blocks.horizon_total_limit[0], generic_sum, rtol=3e-13, atol=3e-13)
            columns = np.arange(1, 8)
            extrapolated = np.tensordot(finite_difference_weights(candidate.R[0], candidate.R[columns], 0),
                                        blocks.total[columns], axes=(0, 0))
            relative = np.linalg.norm(blocks.total[0] - extrapolated) / np.linalg.norm(blocks.total[0])
            assert relative < 0.1


def test_local_source_matches_full_grid_path_over_rk4_steps(plunge):
    grid = build_spatial_grid(0.8, 2, 512, 129)
    config = EvolutionConfig(1, 0.8, 2, 512, 129, flux_model_version=plunge.manifest.flux_model_version)
    trajectory = integrate_trajectory(plunge, config, 4 * grid.dR)
    local_source = GridSourceEvaluator(grid, trajectory, 2, 4 * grid.dR, 4 * grid.dy)
    full_source = GridSourceEvaluator(grid, trajectory, 2, 4 * grid.dR, 4 * grid.dy)
    assert local_source.uses_local_support
    full_source.uses_local_support = False
    for T in (20.0, 100.0, 190.0, trajectory.events.horizon_crossing - 1e-3):
        local, full = local_source.blocks(T), full_source.blocks(T)
        for actual, expected in zip(
            (local.block1, local.block2, local.block3, local.block4, local.total),
            (full.block1, full.block2, full.block3, full.block4, full.total),
        ):
            assert np.allclose(actual, expected, rtol=2e-12, atol=2e-12)
    profile = (1 - grid.yy**2) * np.exp(2j * grid.RR)
    initial = EvolutionState(0.3 * profile, profile)
    dt = 2e-4
    local_final, full_final, T = initial, initial, 100.0
    local_rhs, full_rhs = SpatialRHS(grid, 0.8, 2, local_source), SpatialRHS(grid, 0.8, 2, full_source)
    for _ in range(20):
        local_final = rk4_step(local_rhs, T, local_final, dt)
        full_final = rk4_step(full_rhs, T, full_final, dt)
        T += dt
    for actual, expected in (
        (local_final.P, full_final.P),
        (local_final.psi, full_final.psi),
    ):
        relative = np.linalg.norm(actual - expected) / np.linalg.norm(expected)
        assert relative < 2e-12


def test_rhs_rk4_stage_times_timestep_and_scri_projection(grid, trajectory):
    rhs = SpatialRHS(grid, 0.8, 2, source=None)
    state = zero_state(grid.RR.shape)
    observed = []
    advanced = rk4_step(rhs, 0.0, state, 1e-4, stage_observer=observed.append)
    assert observed == [0.0, 5e-5, 5e-5, 1e-4]
    assert np.array_equal(advanced.P, state.P) and np.array_equal(advanced.psi, state.psi)
    choice = choose_timestep(grid, 0.8, 2, trajectory, 0.1)
    assert choice.dt == pytest.approx(0.1 / choice.steps_per_output)
    assert choice.controlling_limit in {"wave", "motion", "phase", "amplitude", "ramp"}
    sourced_rhs = SpatialRHS(grid, 0.8, 2, GridSourceEvaluator(grid, trajectory, 2, 4 * grid.dR, 4 * grid.dy))
    sourced = rk4_step(sourced_rhs, 10.0, state, choice.dt)
    assert np.all(np.isfinite(sourced.P)) and np.max(np.abs(sourced.P)) > 0
    harmonic = np.sqrt(5 / (64 * np.pi)) * (1 - grid.y) ** 2
    field = np.zeros(grid.RR.shape, dtype=complex); field[-1] = (1.2 - 0.3j) * harmonic
    assert project_scri_mode(field, grid, 2, 2) == pytest.approx(1.2 - 0.3j, rel=2e-5)
    calls = []
    class RecordingSource:
        def __call__(self, T):
            calls.append(T); return np.zeros(grid.RR.shape, complex)
    rk4_step(SpatialRHS(grid, 0.8, 2, RecordingSource()), 1.0, state, 0.2)
    assert calls == [1.0, 1.1, 1.1, 1.2]


def test_timestep_extrema_do_not_underestimate_dense_reference(grid, trajectory):
    choice = choose_timestep(grid, 0.8, 2, trajectory, 0.1, tau_on=20)
    times = np.linspace(choice.t_stable, np.nextafter(trajectory.events.horizon_crossing, -np.inf), 5001)
    rmax = pmax = 0.0
    qmax = 0.0
    pairs = []
    for T in times:
        sample = trajectory.sample(float(T)); rmax = max(rmax, abs(sample.R_jet.first)); pmax = max(pmax, abs(sample.Phi_jet.first))
        phase = np.exp(-2j * sample.Phi)
        for raw in (sample.amplitude_nn, sample.amplitude_mn, sample.amplitude_mm):
            q0 = phase * raw.value
            q1 = phase * (raw.first - 2j * sample.Phi_jet.first * raw.value)
            pairs.append((q0, q1)); qmax = max(qmax, abs(q0))
    gamma = max(abs(q1) / max(abs(q0), choice.q_floor) for q0, q1 in pairs)
    assert 0.5 * (4 * grid.dR) / choice.dt_motion >= rmax
    assert 0.10 / (2 * choice.dt_phase) >= pmax
    assert choice.q_reference >= qmax
    assert choice.gamma_source >= gamma
    assert choice.dt_ramp == pytest.approx(0.5)
    fallback_events = replace(trajectory.events, horizon_crossing=110.0, source_off=110.0)
    fallback = choose_timestep(grid, 0.8, 2, replace(trajectory, events=fallback_events), 0.1, tau_on=20)
    assert fallback.t_stable == 90.0 and fallback.extrema_interval == (90.0, 110.0)
    with pytest.raises(Exception):
        choose_timestep(grid, 0.8, 2,
                        replace(trajectory, events=replace(trajectory.events, horizon_crossing=20.0, source_off=20.0)),
                        0.1, tau_on=20)


def test_checkpoint_restart_mode_buffer_and_memmap(tmp_path):
    def rhs(_T, state):
        return EvolutionState(0.2 * state.P + 1, -0.1 * state.psi + 2)
    initial = EvolutionState(np.zeros((3, 2), complex), np.zeros((3, 2), complex))
    full = evolve_accepted_steps(rhs, initial, 0, 0.4, 0.1)
    half = evolve_accepted_steps(rhs, initial, 0, 0.2, 0.1)
    digest = hash_config("fixed")
    path = tmp_path / "checkpoint.npz"
    save_checkpoint(path, Checkpoint(0.2, half, 0.1, 2, digest, "trajectory-v1", "code-v1"))
    loaded = load_checkpoint(path, expected_shape=(3, 2), expected_config_hash=digest,
                             expected_trajectory_hash="trajectory-v1", expected_code_hash="code-v1")
    resumed = evolve_accepted_steps(rhs, loaded.state, loaded.T, 0.4, loaded.dt)
    assert np.array_equal(resumed.P, full.P) and np.array_equal(resumed.psi, full.psi)
    rolling = RollingCheckpoint(tmp_path / "rolling")
    for index in range(3):
        rolling.save(Checkpoint(0.1 * index, half, 0.1, index, digest, "trajectory-v1", "code-v1"))
    assert len(list((tmp_path / "rolling").glob("checkpoint_*.npz"))) == 2
    assert rolling.latest(expected_shape=(3, 2), expected_config_hash=digest,
                          expected_trajectory_hash="trajectory-v1", expected_code_hash="code-v1").T == pytest.approx(0.2)
    with np.load(path, allow_pickle=False) as data:
        payload = {name: data[name] for name in data.files}
    payload["status"] = "unknown"
    bad_checkpoint = tmp_path / "bad_checkpoint.npz"; np.savez(bad_checkpoint, **payload)
    with pytest.raises(Exception):
        load_checkpoint(bad_checkpoint, expected_shape=(3, 2), expected_config_hash=digest,
                        expected_trajectory_hash="trajectory-v1", expected_code_hash="code-v1")
    for field_name, invalid_count in (("next_output_index", 1.5), ("mode_count", -1),
                                      ("scri_count", "2"), ("field_valid_slices", True)):
        invalid_payload = dict(payload); invalid_payload["status"] = "complete"; invalid_payload[field_name] = invalid_count
        invalid_checkpoint = tmp_path / f"invalid_{field_name}.npz"; np.savez(invalid_checkpoint, **invalid_payload)
        with pytest.raises(Exception):
            load_checkpoint(invalid_checkpoint, expected_shape=(3, 2), expected_config_hash=digest,
                            expected_trajectory_hash="trajectory-v1", expected_code_hash="code-v1")
    mode_metadata = {"ell": 2, "m": 2, "events": {}, "normalization": "test", "error_metadata": {},
                     "input_hash": "input", "code_hash": "code", "config_hash": digest,
                     "analysis_start": 0.0, "phi_dot_at_analysis_start": 0.4,
                     "run_kind": "full", "T_bound": 0.11, "T1": 0.1, "N_end": 1, "delta_t_post": 120.0,
                     "environment": {"python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__}}
    buffer = ModeBuffer(tmp_path / "mode.npz", mode_metadata, 2); buffer.append(0, 1j); buffer.append(0.1, 2j); buffer.flush(complete=True)
    with np.load(tmp_path / "mode.npz") as data:
        assert np.array_equal(data["psi4_lm"], [1j, 2j])
    assert np.array_equal(load_mode(tmp_path / "mode.npz", expected_ell=2, expected_m=2,
                                    expected_environment=mode_metadata["environment"])[1], [1j, 2j])
    with pytest.raises(Exception):
        load_mode(tmp_path / "mode.npz", expected_environment={"python": "tampered"})
    short_buffer = ModeBuffer(tmp_path / "short.npz", mode_metadata, 10)
    short_buffer.append(0, 0j); short_buffer.flush(status="incomplete")
    with pytest.raises(Exception): ModeBuffer.resume(tmp_path / "short.npz", mode_metadata, 10, truncate=2)
    unknown_meta = dict(mode_metadata); unknown_meta.update({"schema": "c1-mode-v1", "status": "unknown"})
    with (tmp_path / "unknown_mode.npz").open("wb") as stream:
        np.savez_compressed(stream, T=np.array([0.0]), psi4_lm=np.array([0j]), metadata=json.dumps(unknown_meta))
    with pytest.raises(Exception): load_mode(tmp_path / "unknown_mode.npz", allow_incomplete=True)
    _, _, valid_mode_meta = load_mode(tmp_path / "mode.npz", allow_incomplete=True)
    for invalid_count in (1.5, -1, "2", True):
        invalid_meta = dict(valid_mode_meta); invalid_meta["sample_count"] = invalid_count
        invalid_path = tmp_path / f"invalid_mode_{str(invalid_count).replace('.', '_')}.npz"
        with invalid_path.open("wb") as stream:
            np.savez_compressed(stream, T=np.array([0.0, 0.1]), psi4_lm=np.array([1j, 2j]),
                                metadata=json.dumps(invalid_meta))
        with pytest.raises(Exception): load_mode(invalid_path, allow_incomplete=True)
        invalid_base = dict(mode_metadata); invalid_base["N_end"] = invalid_count
        with pytest.raises(Exception): ModeBuffer(tmp_path / "never_written.npz", invalid_base, 2)
    buffer.append(0.2, 3j); buffer.flush(status="incomplete")
    resumed_buffer = ModeBuffer.resume(tmp_path / "mode.npz", mode_metadata, 2, truncate=2)
    assert len(resumed_buffer.times) == 2
    assert len(load_mode(tmp_path / "mode.npz", allow_incomplete=True)[0]) == 2
    memmap = FullFieldMemmap(tmp_path / "field.npy", 2, (3, 2)); memmap.append(full.psi); memmap.flush()
    assert memmap.valid_slices == 1 and np.array_equal(np.load(tmp_path / "field.npy")[0], full.psi)


def test_config_binding_driver_restart_and_outputs(tmp_path, plunge):
    cfg = EvolutionConfig(1, 0.8, 2, 128, 33, sigma_r_over_dr=3.5, sigma_y_over_dy=3.0,
                          flux_model_version=plunge.manifest.flux_model_version, tau_on_over_m=8)
    output = OutputRequest((2,), save_field=True, save_scri_field=True, output_dt_over_m=5e-5,
                           field_dt_over_m=1e-4, checkpoint_dt_over_m=1e-4, chunk_size=2)
    system = assemble_single_system(cfg, plunge, output)
    assert system.source.sigma_R == pytest.approx(3.5 * system.grid.dR)
    assert system.source.sigma_y == pytest.approx(3.0 * system.grid.dy)
    assert system.source.tau_on == 8 and system.output == output
    assert system.N_end == math.floor(system.T_bound / output.output_dt_over_m)
    assert system.T1 == pytest.approx(system.N_end * output.output_dt_over_m)
    assert system.T1 <= system.T_bound and len(system.code_hash) == 64
    packages = (
        Path(__file__).parents[1] / "src" / "sminus2_point_particle",
        Path(kerr_waveform_tools.__file__).resolve().parent,
    )
    digest = hashlib.sha256()
    for package in packages:
        for path in sorted(package.glob("*.py")):
            digest.update(package.name.encode("ascii"))
            digest.update(path.name.encode("ascii"))
            digest.update(path.read_bytes())
    assert system.code_hash == digest.hexdigest()
    with pytest.raises(Exception):
        assemble_single_system(replace(cfg, chi=0.7), plunge, output)
    run_dir = tmp_path / "run"
    first = run_single_system(system, run_dir, T_end=2.25e-4)
    assert first.final_time == pytest.approx(2.25e-4)
    first_T = load_mode(run_dir / "psi4_l2_m2.npz", expected_ell=2, expected_m=2, allow_incomplete=True)[0]
    assert np.allclose(np.diff(first_T), output.output_dt_over_m)
    assert first_T[-1] == pytest.approx(2e-4)
    resumed = run_single_system(system, run_dir, T_end=3e-4, restart=True)
    assert resumed.final_time == pytest.approx(3e-4)
    T, mode, metadata = load_mode(run_dir / "psi4_l2_m2.npz", expected_ell=2, expected_m=2, allow_incomplete=True)
    assert metadata["status"] == "test" and len(T) == resumed.output_count
    assert np.allclose(np.diff(T), output.output_dt_over_m) and T[-1] == pytest.approx(3e-4)
    assert (run_dir / "scri_field.npz").exists()
    run_metadata = json.loads((run_dir / "metadata.json").read_text(encoding="ascii"))
    assert run_metadata["status"] == "test" and run_metadata["output"]["chunk_size"] == 2
    assert (run_metadata["T_bound"], run_metadata["T1"], run_metadata["N_end"], run_metadata["delta_t_post"]) == pytest.approx(
        (system.T_bound, system.T1, system.N_end, 120.0))
    assert run_metadata["valid_field_slices"] == 3
    assert np.load(run_dir / "field.npy", mmap_mode="r").shape[0] >= 3
    assert run_metadata["environment"] == {"python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__}
    scri_T, scri_field, scri_metadata = load_scri_field(
        run_dir / "scri_field.npz", expected_config_hash=system.config_hash,
        expected_trajectory_hash=system.trajectory_hash, expected_code_hash=system.code_hash,
        allow_incomplete=True)
    assert len(scri_T) == 6 and scri_field.shape == (6, 33) and scri_metadata["status"] == "test"
    with pytest.raises(Exception): load_scri_field(run_dir / "scri_field.npz")
    with pytest.raises(Exception): load_scri_field(run_dir / "scri_field.npz", expected_code_hash="tampered", allow_incomplete=True)
    incomplete_scri = run_dir / "incomplete_scri.npz"
    incomplete_meta = dict(scri_metadata); incomplete_meta["status"] = "incomplete"
    with incomplete_scri.open("wb") as stream:
        np.savez_compressed(stream, T=scri_T, psi4_m=scri_field, metadata=json.dumps(incomplete_meta))
    with pytest.raises(Exception): load_scri_field(incomplete_scri)
    unknown_scri_meta = dict(scri_metadata); unknown_scri_meta["status"] = "unknown"
    unknown_scri = run_dir / "unknown_scri.npz"
    with unknown_scri.open("wb") as stream:
        np.savez_compressed(stream, T=scri_T, psi4_m=scri_field, metadata=json.dumps(unknown_scri_meta))
    with pytest.raises(Exception): load_scri_field(unknown_scri, allow_incomplete=True)
    for invalid_count in (1.5, -1, "6", True):
        invalid_meta = dict(scri_metadata); invalid_meta["sample_count"] = invalid_count
        invalid_path = run_dir / f"invalid_scri_{str(invalid_count).replace('.', '_')}.npz"
        with invalid_path.open("wb") as stream:
            np.savez_compressed(stream, T=scri_T, psi4_m=scri_field, metadata=json.dumps(invalid_meta))
        with pytest.raises(Exception): load_scri_field(invalid_path, allow_incomplete=True)
    full_field, field_metadata = load_full_field(
        run_dir / "field.npy", run_dir / "field_metadata.json", expected_config_hash=system.config_hash,
        expected_trajectory_hash=system.trajectory_hash, expected_code_hash=system.code_hash,
        allow_incomplete=True)
    assert full_field.shape == (3,) + system.grid.RR.shape and field_metadata["valid_slices"] == 3
    with pytest.raises(Exception): load_full_field(run_dir / "field.npy", run_dir / "field_metadata.json")
    bad_metadata = dict(field_metadata); bad_metadata["valid_slices"] = 99
    bad_path = run_dir / "bad_field_metadata.json"; bad_path.write_text(json.dumps(bad_metadata), encoding="ascii")
    with pytest.raises(Exception): load_full_field(run_dir / "field.npy", bad_path, allow_incomplete=True)
    incomplete_field = dict(field_metadata); incomplete_field["status"] = "incomplete"
    incomplete_path = run_dir / "incomplete_field_metadata.json"; incomplete_path.write_text(json.dumps(incomplete_field), encoding="ascii")
    with pytest.raises(Exception): load_full_field(run_dir / "field.npy", incomplete_path)
    unknown_field = dict(field_metadata); unknown_field["status"] = "unknown"
    unknown_field_path = run_dir / "unknown_field_metadata.json"; unknown_field_path.write_text(json.dumps(unknown_field), encoding="ascii")
    with pytest.raises(Exception): load_full_field(run_dir / "field.npy", unknown_field_path, allow_incomplete=True)
    for invalid_count in (1.5, -1, "3", True):
        invalid_meta = dict(field_metadata); invalid_meta["valid_slices"] = invalid_count
        invalid_path = run_dir / f"invalid_field_{str(invalid_count).replace('.', '_')}.json"
        invalid_path.write_text(json.dumps(invalid_meta), encoding="ascii")
        with pytest.raises(Exception): load_full_field(run_dir / "field.npy", invalid_path, allow_incomplete=True)
    assert load_run_metadata(run_dir, allow_incomplete=True)["T_end"] == pytest.approx(3e-4)
    checkpoint = RollingCheckpoint(run_dir / "checkpoints").latest(
        expected_shape=system.grid.RR.shape, expected_config_hash=system.config_hash,
        expected_trajectory_hash=system.trajectory_hash, expected_code_hash=system.code_hash)
    assert (checkpoint.T, checkpoint.mode_count, checkpoint.scri_count, checkpoint.field_valid_slices) == pytest.approx((3e-4, 6, 6, 3))
    mode_short = tmp_path / "mode-short"; shutil.copytree(run_dir, mode_short)
    mode_T, mode_values, mode_meta = load_mode(mode_short / "psi4_l2_m2.npz", allow_incomplete=True)
    with (mode_short / "psi4_l2_m2.npz").open("wb") as stream:
        np.savez_compressed(stream, T=mode_T[:-1], psi4_lm=mode_values[:-1], metadata=json.dumps(mode_meta))
    with pytest.raises(Exception): run_single_system(system, mode_short, T_end=3.5e-4, restart=True)
    scri_short = tmp_path / "scri-short"; shutil.copytree(run_dir, scri_short)
    with (scri_short / "scri_field.npz").open("wb") as stream:
        np.savez_compressed(stream, T=scri_T[:-1], psi4_m=scri_field[:-1], metadata=json.dumps(scri_metadata))
    with pytest.raises(Exception): run_single_system(system, scri_short, T_end=3.5e-4, restart=True)
    field_short = tmp_path / "field-short"; shutil.copytree(run_dir, field_short)
    short_field_metadata = dict(field_metadata); short_field_metadata["valid_slices"] = 2
    (field_short / "field_metadata.json").write_text(json.dumps(short_field_metadata), encoding="ascii")
    with pytest.raises(Exception): run_single_system(system, field_short, T_end=3.5e-4, restart=True)
    tampered_run = dict(run_metadata); tampered_run["environment"] = {"python": "tampered"}
    (run_dir / "metadata.json").write_text(json.dumps(tampered_run), encoding="ascii")
    with pytest.raises(Exception): run_single_system(system, run_dir, T_end=3.5e-4, restart=True)
    tampered_run["status"] = "unknown"
    (run_dir / "metadata.json").write_text(json.dumps(tampered_run), encoding="ascii")
    with pytest.raises(Exception): load_run_metadata(run_dir, allow_incomplete=True)
    for field_name, invalid_count in (("N_end", 1.5), ("output_count", -1), ("valid_field_slices", True)):
        invalid_run = dict(run_metadata); invalid_run[field_name] = invalid_count
        (run_dir / "metadata.json").write_text(json.dumps(invalid_run), encoding="ascii")
        with pytest.raises(Exception): load_run_metadata(run_dir, allow_incomplete=True)


def test_real_short_driver_mode_uses_absolute_science_subgrid(tmp_path, plunge):
    config = EvolutionConfig(1, 0.8, 2, 128, 33, flux_model_version=plunge.manifest.flux_model_version)
    output = OutputRequest((2,), output_dt_over_m=0.05, field_dt_over_m=0.1,
                           checkpoint_dt_over_m=0.1, chunk_size=8)
    system = assemble_single_system(config, plunge, output)
    run_single_system(system, tmp_path / "absolute-grid", T_end=0.25)
    T, _, metadata = load_mode(tmp_path / "absolute-grid" / "psi4_l2_m2.npz",
                               expected_ell=2, expected_m=2, allow_incomplete=True)
    assert metadata["status"] == "test" and T[0] == pytest.approx(0.05)
    indices = absolute_coarse_indices(T)
    assert np.allclose(T[indices], [0.1, 0.2])
    comparison_T = np.arange(T[0], 120.0001, 0.05)
    omega = 0.8
    comparison = compare_sampling(comparison_T, -omega**2 * np.exp(1j * omega * comparison_T),
                                  ell=2, m=2, phi_dot_at_start=0.4, epsilon_v2=0.02,
                                  analysis_start=10.1, stop=120.0)
    assert comparison.passed


def test_mode_to_convert_strain_cli_schema_and_provenance(tmp_path):
    T = np.arange(0, 120.0001, 0.05); omega = 0.8
    psi4 = -omega**2 * np.exp(1j * omega * T)
    metadata = {"ell": 2, "m": 2, "events": {"source_off": 1.0}, "normalization": "peeling psi4",
                "error_metadata": {"epsilon_v2": 0.02}, "input_hash": "input-hash",
                "code_hash": "code-hash", "config_hash": "config-hash", "analysis_start": 100.0,
                "phi_dot_at_analysis_start": 0.4, "run_kind": "test", "T_bound": 120.03,
                "T1": 120.0, "N_end": 2400, "delta_t_post": 120.0,
                "environment": {"python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__}}
    metadata["background"] = {"mass_scale": 1.0, "M_internal": 1.0, "chi": 0.8,
                              "a_internal": 0.8, "L_internal": 1.0}
    metadata["amplitude_scaling"] = "unit particle-mass response"
    mode_path, output_path = tmp_path / "mode.npz", tmp_path / "H.npz"
    buffer = ModeBuffer(mode_path, metadata, 10000)
    for time, value in zip(T, psi4): buffer.append(time, value)
    buffer.flush(status="test")
    rejected = subprocess.run([sys.executable, "scripts/convert_strain.py", str(mode_path), str(output_path),
                               "--ell", "2", "--m", "2"], capture_output=True)
    assert rejected.returncode != 0
    subprocess.run([sys.executable, "scripts/convert_strain.py", str(mode_path), str(output_path),
                    "--ell", "2", "--m", "2", "--phi-dot-at-start", "0.4", "--test-input",
                    "--cutoff-factor", "0.6", "--taper-width", "9"], check=True)
    with np.load(output_path, allow_pickle=False) as data:
        converted = json.loads(str(data["metadata"]))
        assert data["H_lm"].shape == T.shape
    assert converted["schema"] == "c2-H-v2" and converted["status"] == "test"
    assert converted["output_object"] == "H_lm=r*h_lm"
    assert converted["source_mode_metadata"]["input_hash"] == "input-hash"
    assert converted["analysis_start"] == 100 and converted["analysis_stop"] == 120
    assert converted["transform_interval"] == [T[0], 120]
    assert converted["analysis_interval"] == [100, 120]
    assert converted["recommended_untapered_interval"] == [100, 111]
    assert converted["cutoff_factor"] == 0.6 and converted["taper_width"] == 9
    assert converted["input_schema"] == "c1-mode-v2" and converted["sampling_status"] == "owner-waived"
    assert np.array_equal(np.load(output_path)["taper"], quintic_taper(T, T[0], 120, 9))
    assert np.load(output_path)["taper"][np.searchsorted(T, 50.0)] == 1.0
    assert "forward exp" in converted["fft_convention"]
    assert converted["source_environment"] == metadata["environment"]
    loaded, loaded_metadata = load_ffi_result(output_path, expected_ell=2, expected_m=2, allow_test=True)
    assert np.array_equal(loaded.T, T) and loaded_metadata["status"] == "test"


def test_ffi_bridge_taper_and_sampling_candidate():
    T = np.arange(0, 120.0001, 0.05)
    omega = 0.8
    H = np.exp(1j * omega * T)
    psi4 = -omega**2 * H
    result = fixed_frequency_integrate(T, psi4, ell=2, m=2, phi_dot_at_start=omega / 2, analysis_start=10, stop=120)
    assert result.metadata.output_object == "H_lm=r*h_lm"
    assert result.metadata.L2_bridge == 1 and result.metadata.omega0 == pytest.approx(0.75 * omega)
    assert result.metadata.a1_status == "conditional-open"
    assert len(result.metadata.input_hash) == 64
    window = quintic_taper(T, T[0], 120, 10)
    assert window[0] == 0 and np.max(window) == 1
    interior = (T >= 30) & (T <= 100)
    phase_aligned = np.vdot(H[interior], result.H[interior]) / np.vdot(H[interior], H[interior])
    assert phase_aligned == pytest.approx(1, rel=0.03, abs=0.03)
    comparison = compare_sampling(T, psi4, ell=2, m=2, phi_dot_at_start=omega / 2, epsilon_v2=0.02, analysis_start=10, stop=120)
    assert comparison.threshold == 1e-3
    assert comparison.passed and comparison.epsilon < comparison.threshold
    low_cutoff = fixed_frequency_integrate(T, psi4, ell=2, m=2, phi_dot_at_start=omega / 2,
                                           analysis_start=10, stop=120, cutoff_factor=0.5)
    high_cutoff = fixed_frequency_integrate(T, psi4, ell=2, m=2, phi_dot_at_start=omega / 2,
                                            analysis_start=10, stop=120, cutoff_factor=1.0)
    assert np.all(np.isfinite(low_cutoff.H)) and np.all(np.isfinite(high_cutoff.H))


@pytest.mark.parametrize(
    "save_psi4,save_strain",
    [(False, False), (True, False), (False, True), (True, True)],
)
def test_finalize_mode_products_honors_output_request(
    tmp_path,
    plunge,
    save_psi4,
    save_strain,
):
    output = OutputRequest(
        (2,),
        save_psi4_lm=save_psi4,
        save_strain_lm=save_strain,
        output_dt_over_m=0.05,
    )
    config = EvolutionConfig(
        1,
        0.8,
        2,
        128,
        33,
        flux_model_version=plunge.manifest.flux_model_version,
    )
    system = replace(
        assemble_single_system(config, plunge, output),
        T1=120.0,
        N_end=2400,
    )
    directory = tmp_path / f"{int(save_psi4)}-{int(save_strain)}"
    directory.mkdir()
    mode_path = directory / "psi4_l2_m2.npz"
    if save_psi4 or save_strain:
        T = np.arange(0.0, 120.0001, 0.05)
        omega = 0.8
        metadata = {
            "ell": 2,
            "m": 2,
            "events": asdict(system.trajectory.events),
            "normalization": "test peeling psi4",
            "error_metadata": {},
            "input_hash": "input-hash",
            "code_hash": system.code_hash,
            "config_hash": system.config_hash,
            "analysis_start": 100.0,
            "phi_dot_at_analysis_start": omega / 2,
            "run_kind": "full",
            "T_bound": 120.0,
            "T1": 120.0,
            "N_end": 2400,
            "delta_t_post": 120.0,
            "environment": {"python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__},
            "background": {"mass_scale": 1.0, "M_internal": 1.0, "chi": 0.8,
                           "a_internal": 0.8, "L_internal": 1.0},
            "amplitude_scaling": "unit particle-mass response",
        }
        buffer = ModeBuffer(mode_path, metadata, len(T) + 1)
        for time, value in zip(T, -omega**2 * np.exp(1j * omega * T)):
            buffer.append(time, value)
        buffer.flush(status="complete")
    finalize_mode_products(system, directory, full_run=True)
    assert mode_path.exists() is save_psi4
    assert (directory / "psi4_l2_m2_real_abs.png").exists() is save_psi4
    assert (directory / "H_l2_m2.npz").exists() is save_strain
    assert (directory / "H_l2_m2_real_abs.png").exists() is save_strain
    if save_strain:
        load_ffi_result(directory / "H_l2_m2.npz", expected_ell=2, expected_m=2)


def test_short_run_rejects_requested_strain_before_writing(tmp_path, plunge):
    output = OutputRequest((2,), save_psi4_lm=False, save_strain_lm=True)
    config = EvolutionConfig(
        1,
        0.8,
        2,
        128,
        33,
        flux_model_version=plunge.manifest.flux_model_version,
    )
    system = assemble_single_system(config, plunge, output)
    directory = tmp_path / "short-strain"
    with pytest.raises(ValueError, match="complete run"):
        run_single_system(system, directory, T_end=1e-4)
    assert not directory.exists()


def test_run_case_consumes_initial_data_cache_miss_then_hit(tmp_path, plunge):
    config = EvolutionConfig(
        1,
        0.8,
        2,
        128,
        33,
        flux_model_version=plunge.manifest.flux_model_version,
    )
    output = OutputRequest((2,), output_dt_over_m=5e-5, checkpoint_dt_over_m=1e-4,
                           field_dt_over_m=1e-4, chunk_size=8)
    config_path, output_path = tmp_path / "config.json", tmp_path / "output.json"
    config_path.write_text(config.to_json(), encoding="ascii")
    output_path.write_text(output.to_json(), encoding="ascii")
    cache = tmp_path / "initial-data"
    statuses = []
    keys = []
    for index in range(2):
        directory = tmp_path / f"run-{index}"
        subprocess.run(
            [
                sys.executable,
                "scripts/run_case.py",
                str(config_path),
                str(output_path),
                str(FIXTURE),
                str(directory),
                "--initial-data-cache",
                str(cache),
                "--t-end",
                "0.0001",
                "--no-progress",
            ],
            check=True,
        )
        metadata = load_run_metadata(directory, allow_incomplete=True)
        statuses.append(metadata["initial_data_provenance"]["status"])
        keys.append(metadata["initial_data_provenance"]["key"])
    assert statuses == ["miss-generated", "hit"]
    assert keys[0] == keys[1]
    assert len(list(cache.glob("*.json"))) == 1


def test_direct_data_staging_finalize_and_missing_product_guard(tmp_path):
    config = EvolutionConfig(1, 0.8, 2, 128, 33, flux_model_version="storage-test")
    no_products = OutputRequest((2,), save_psi4_lm=False, save_strain_lm=False)
    paths = data_run_paths(tmp_path / "runs", "complete-case")
    (paths.staging / "checkpoints").mkdir(parents=True)
    (paths.staging / "checkpoints" / "temporary").write_text("restart", encoding="ascii")
    metadata = {
        "schema": "c1-run-v2",
        "run_kind": "full",
        "config_hash": "config",
        "code_hash": "code",
        "trajectory_hash": "trajectory",
    }
    finalize_run(paths.staging, metadata, status="complete")
    final = finalize_data_bundle(paths, config, no_products)
    assert final == paths.final and final.is_dir() and not paths.staging.exists()
    assert not (final / "checkpoints").exists()
    bundle = json.loads((final / "bundle.json").read_text(encoding="ascii"))
    assert bundle["status"] == "complete" and bundle["review_status"] == "awaiting_review"

    guarded = data_run_paths(tmp_path / "runs", "missing-products")
    finalize_run(guarded.staging, metadata, status="complete")
    with pytest.raises(Exception, match="mode"):
        finalize_data_bundle(guarded, config, OutputRequest((2,), save_psi4_lm=True))
    assert guarded.staging.exists() and not guarded.final.exists()


def test_output_defaults_match_plan():
    output = OutputRequest()
    assert (output.output_dt_over_m, output.field_dt_over_m, output.checkpoint_dt_over_m) == (0.1, 5.0, 10.0)
    with pytest.raises(Exception): EvolutionConfig(1, 0.8, 3, 128, 33, flux_model_version="v")
    assert EvolutionConfig(1, 0.8, 2, 128, 33, flux_model_version="v", tau_on_over_m=21).tau_on_over_m == 21
    with pytest.raises(Exception): EvolutionConfig(1, 0.8, 2, 128, 33, flux_model_version="v", tau_on_over_m=41)
    with pytest.raises(Exception): quintic_window(1, 41)
    with pytest.raises(Exception): GridSourceEvaluator(grid, trajectory, 2, 4 * grid.dR, 4 * grid.dy, tau_on=41)
    with pytest.raises(Exception): choose_timestep(grid, 0.8, 2, trajectory, tau_on=-1)
    with pytest.raises(Exception): choose_timestep(grid, 0.8, 2, trajectory, tau_on=41)


def test_m4_default_output_and_explicit_invalid_ell():
    flux = FluxRecord(0.8, 6, 0.01687986787530486, 0.01, "m4-fixture", "synthetic m4 interface fixture")
    initial = generate_ori_thorne_initial_data(flux)
    config = EvolutionConfig(1, 0.8, 4, 128, 33, flux_model_version=flux.flux_version)
    system = assemble_single_system(config, initial)
    assert system.output.ell_out == (4,)
    field = np.zeros(system.grid.RR.shape, complex)
    assert np.isfinite(project_scri_mode(field, system.grid, 4, 4))
    with pytest.raises(Exception):
        assemble_single_system(config, initial, OutputRequest(ell_out=(2,)))
