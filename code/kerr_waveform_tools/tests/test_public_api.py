import json

import numpy as np
import pytest

from kerr_waveform_tools import (
    ArtifactStatus,
    FluxRecord,
    JsonInitialDataCache,
    ModeBuffer,
    TransitionSpec,
    fixed_frequency_integrate,
    generate_ori_thorne_initial_data,
    integrate_trajectory,
    isco_quantities,
    load_mode,
    project_mode,
    projection_row,
    quintic_source_turn_on,
    quintic_source_turn_on_values,
)


def test_isco_transition_key_and_json_cache(tmp_path):
    flux = FluxRecord(
        0.8,
        6,
        0.01687986787530486,
        0.01660590066635954,
        "public-test-l6",
        "recorded regression fixture",
    )
    initial = generate_ori_thorne_initial_data(flux)
    spec = TransitionSpec(0.8, flux.flux_version)
    key = spec.initial_data_key(
        initial.manifest.schema_version,
        initial.manifest.generator_version,
    )
    cache = JsonInitialDataCache(tmp_path, spec)
    assert cache.load(key) is None
    cache.store(key, initial)
    assert cache.load(key) == initial
    assert isco_quantities(0.8)[0] == pytest.approx(2.9066438544641957)


def test_worldline_events_are_ordered():
    flux = FluxRecord(
        0.8,
        6,
        0.01687986787530486,
        0.01660590066635954,
        "public-worldline-l6",
        "recorded regression fixture",
    )
    initial = generate_ori_thorne_initial_data(flux)
    trajectory = integrate_trajectory(
        initial,
        TransitionSpec(0.8, flux.flux_version),
        sigma_R=0.02,
    )
    assert (
        trajectory.events.light_ring
        < trajectory.events.horizon_crossing
        == trajectory.events.source_off
        < trajectory.events.terminal
    )
    jet = trajectory.sample(100.0).R_jet
    assert np.all(np.isfinite([jet.value, jet.first, jet.second]))


def test_spherical_projection_recovers_mode():
    y = np.linspace(-1, 1, 1001)
    weights = np.full_like(y, 2 / (len(y) - 1))
    weights[[0, -1]] *= 0.5
    row = projection_row(2, 2, y, weights)
    harmonic = np.sqrt(5 / (64 * np.pi)) * (1 - y) ** 2
    assert project_mode((1.2 - 0.3j) * harmonic, row) == pytest.approx(
        1.2 - 0.3j,
        rel=2e-6,
    )


def test_shared_source_turn_on_value_and_jets():
    off = quintic_source_turn_on(0.0)
    on = quintic_source_turn_on(20.0)
    assert (off.value, off.first, off.second) == (0.0, 0.0, 0.0)
    assert (on.value, on.first, on.second) == (1.0, 0.0, 0.0)
    middle = quintic_source_turn_on(10.0)
    assert middle.value == pytest.approx(0.5)
    np.testing.assert_allclose(
        quintic_source_turn_on_values(np.array([0.0, 10.0, 20.0])),
        [0.0, 0.5, 1.0],
    )


def test_mode_io_and_full_series_ffi(tmp_path):
    T = np.arange(0.0, 120.0001, 0.05)
    omega = 0.8
    psi4 = -omega**2 * np.exp(1j * omega * T)
    metadata = {
        "ell": 2,
        "m": 2,
        "events": {"source_off": 80.0},
        "normalization": "test",
        "error_metadata": {},
        "input_hash": "input",
        "code_hash": "code",
        "config_hash": "config",
        "analysis_start": 100.0,
        "phi_dot_at_analysis_start": 0.4,
        "run_kind": "test",
        "T_bound": 120.0,
        "T1": 120.0,
        "N_end": len(T) - 1,
        "delta_t_post": 120.0,
        "environment": {"python": "test"},
    }
    path = tmp_path / "mode.npz"
    buffer = ModeBuffer(path, metadata, chunk_size=len(T) + 1)
    for time, value in zip(T, psi4):
        buffer.append(time, value)
    buffer.flush(status=ArtifactStatus.TEST.value)
    stored_T, stored, stored_metadata = load_mode(path, allow_incomplete=True)
    assert np.array_equal(stored_T, T)
    assert np.array_equal(stored, psi4)
    assert stored_metadata["status"] == "test"

    result = fixed_frequency_integrate(
        stored_T,
        stored,
        ell=2,
        m=2,
        phi_dot_at_start=0.4,
        analysis_start=100.0,
        stop=120.0,
        taper_width=10.0,
    )
    assert result.metadata.transform_interval == (0.0, 120.0)
    assert result.metadata.analysis_interval == (100.0, 120.0)
    assert np.all(np.isfinite(result.H))
