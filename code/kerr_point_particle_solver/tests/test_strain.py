import json
import hashlib
from pathlib import Path

import numpy as np
import pytest

from sminus2_point_particle.errors import ContractError
from sminus2_point_particle.strain import (
    compare_cutoff_variation,
    fixed_frequency_integrate,
    load_ffi_result,
    quintic_taper,
    write_ffi_result,
)


def analytic_series():
    T = np.arange(0.0, 200.0001, 0.05)
    omega = 0.8
    H = np.exp(1j * omega * T)
    return T, -omega**2 * H, H, omega


def test_ffi_preserves_full_series_taper_and_metadata():
    T, psi4, H, omega = analytic_series()
    result = fixed_frequency_integrate(T, psi4, ell=2, m=2, phi_dot_at_start=omega / 2,
                                       analysis_start=100, stop=200, L=1.3)
    assert np.array_equal(result.T, T)
    assert result.H.shape == result.taper.shape == T.shape
    assert np.array_equal(result.taper, quintic_taper(T, T[0], 200, 10))
    assert result.metadata.analysis_start == 100 and result.metadata.analysis_stop == 200
    assert (result.metadata.recommended_start, result.metadata.recommended_stop) == (100, 190)
    assert result.metadata.transform_interval == (T[0], 200)
    assert result.metadata.analysis_interval == (100, 200)
    assert result.metadata.L2_bridge == pytest.approx(1.3**2)
    assert result.metadata.sampling_status == "owner-waived"
    assert result.metadata.a1_status == "conditional-open"
    assert len(result.metadata.implementation_hash) == 64
    assert set(result.metadata.implementation_environment) == {"python", "numpy"}
    assert "forward exp" in result.metadata.fft_convention


@pytest.mark.parametrize("mutation", [
    "nan_time", "decreasing", "duplicate", "nonuniform", "nan_real", "inf_imag",
    "start_below", "stop_above", "reversed_interval", "zero_L", "zero_cutoff",
    "zero_taper", "wide_taper", "nan_phase", "zero_frequency", "bad_mode", "wrong_transform_start",
    "analysis_in_left_taper",
])
def test_ffi_rejects_invalid_inputs(mutation):
    T, psi4, _, omega = analytic_series()
    T, psi4 = T.copy(), psi4.copy()
    kwargs = dict(ell=2, m=2, phi_dot_at_start=omega / 2, analysis_start=100, stop=200)
    if mutation == "nan_time": T[4] = np.nan
    elif mutation == "decreasing": T[4] = T[3] - 0.1
    elif mutation == "duplicate": T[4] = T[3]
    elif mutation == "nonuniform": T[4] += 1e-4
    elif mutation == "nan_real": psi4[4] = np.nan + 1j
    elif mutation == "inf_imag": psi4[4] = 1 + 1j * np.inf
    elif mutation == "start_below": kwargs["analysis_start"] = -1
    elif mutation == "stop_above": kwargs["stop"] = 201
    elif mutation == "reversed_interval": kwargs.update(analysis_start=180, stop=170)
    elif mutation == "zero_L": kwargs["L"] = 0
    elif mutation == "zero_cutoff": kwargs["cutoff_factor"] = 0
    elif mutation == "zero_taper": kwargs["taper_width"] = 0
    elif mutation == "wide_taper": kwargs["taper_width"] = 101
    elif mutation == "nan_phase": kwargs["phi_dot_at_start"] = np.nan
    elif mutation == "zero_frequency": kwargs["phi_dot_at_start"] = 0
    elif mutation == "bad_mode": kwargs.update(ell=2, m=4)
    elif mutation == "wrong_transform_start": kwargs["transform_start"] = T[1]
    elif mutation == "analysis_in_left_taper": kwargs["analysis_start"] = 5
    with pytest.raises(ContractError):
        fixed_frequency_integrate(T, psi4, **kwargs)


def test_pre_analysis_signal_participates_in_full_transform():
    T = np.arange(0.05, 200.0001, 0.05)
    psi4 = np.zeros_like(T, dtype=complex)
    early = (T >= 20) & (T <= 40)
    psi4[early] = np.exp(0.7j * T[early])
    result = fixed_frequency_integrate(T, psi4, ell=2, m=2, phi_dot_at_start=0.4,
                                       analysis_start=100, stop=T[-1])
    assert result.taper[np.argmin(np.abs(T - 20))] == pytest.approx(1.0)
    assert np.linalg.norm(result.H) > 0
    assert result.metadata.transform_start == pytest.approx(0.05)
    assert result.metadata.analysis_start == 100


def test_retained_frequency_forward_reconstruction():
    T, psi4, _, omega0 = analytic_series()
    result = fixed_frequency_integrate(T, psi4, ell=2, m=2, phi_dot_at_start=omega0 / 2,
                                       analysis_start=100, stop=200)
    frequencies = 2 * np.pi * np.fft.fftfreq(len(T), d=T[1] - T[0])
    retained = np.abs(frequencies) >= result.metadata.omega0
    reconstructed = -(frequencies[retained] ** 2) * np.fft.fft(result.H)[retained]
    target = np.fft.fft(result.taper * psi4)[retained]
    assert np.linalg.norm(reconstructed - target) / np.linalg.norm(target) < 3e-12


def test_p19_style_low_frequency_pollution_is_cutoff_bounded():
    T, psi4, _, omega = analytic_series()
    low_omega = 0.03
    contamination = 2e-4 * np.exp(1j * low_omega * T)
    polluted = psi4 + contamination
    result = fixed_frequency_integrate(T, polluted, ell=2, m=2, phi_dot_at_start=omega / 2,
                                       analysis_start=100, stop=200)
    assert np.all(np.isfinite(result.H))
    # FFI bounds the nominal low-frequency response by A/omega0^2 rather than A/low_omega^2.
    assert np.max(np.abs(result.H)) < 2.0
    assert 2e-4 / result.metadata.omega0**2 < 2e-4 / low_omega**2


def test_cutoff_variation_is_callable_diagnostic():
    T, psi4, _, omega = analytic_series()
    diagnostic = compare_cutoff_variation(T, psi4, ell=2, m=2, phi_dot_at_start=omega / 2,
                                          analysis_start=100, stop=200)
    assert diagnostic.factors == (0.5, 0.75, 1.0)
    assert diagnostic.relative_differences[1] == pytest.approx(0.0)
    assert np.all(np.isfinite(diagnostic.relative_differences))


def test_atomic_writer_and_loader_reject_tamper(tmp_path):
    T, psi4, _, omega = analytic_series()
    background = {"mass_scale": 1.0, "M_internal": 1.0, "chi": 0.8,
                  "a_internal": 0.8, "L_internal": 1.0}
    result = fixed_frequency_integrate(T, psi4, ell=2, m=2, phi_dot_at_start=omega / 2,
                                       analysis_start=100, stop=200, source_input_hash="source-hash",
                                       source_events={"source_off": 10}, source_normalization="peeling",
                                       source_error_metadata={"V2": "not_run"}, source_environment={"python": "test"},
                                       source_background=background, amplitude_scaling="unit response",
                                       input_schema="c1-mode-v1")
    path = tmp_path / "H.npz"
    source_metadata = {"schema": "c1-mode-v1", "status": "test", "ell": 2, "m": 2,
                       "events": {"source_off": 10}, "normalization": "peeling",
                       "error_metadata": {"V2": "not_run"}, "input_hash": "source-hash",
                       "environment": {"python": "test"}, "run_kind": "test", "T1": 200.0,
                       "analysis_start": 100.0, "phi_dot_at_analysis_start": omega / 2,
                       "background": background, "amplitude_scaling": "unit response"}
    write_ffi_result(path, result, status="test", source_mode_metadata=source_metadata)
    assert path.exists() and not path.with_suffix(".npz.tmp").exists()
    with pytest.raises(ContractError): load_ffi_result(path)
    loaded, metadata = load_ffi_result(path, expected_ell=2, expected_m=2, allow_test=True)
    assert np.array_equal(loaded.T, T) and np.array_equal(loaded.taper, result.taper)
    assert metadata["sample_count"] == len(T) and metadata["sampling_status"] == "owner-waived"

    with np.load(path, allow_pickle=False) as data:
        arrays = {name: data[name] for name in data.files}
    original_metadata = json.loads(str(arrays["metadata"]))
    for name, value in (("schema", "unknown"), ("status", "unknown"),
                        ("sample_count", 1.5), ("output_hash", "tampered"),
                        ("ell", 3), ("analysis_stop", 199.0), ("omega0", 0.123),
                        ("sampling_epsilon", 1e-4), ("implementation_hash", "unknown")):
        changed = dict(original_metadata); changed[name] = value
        tampered = tmp_path / f"tampered_{name}.npz"
        with tampered.open("wb") as stream:
            np.savez_compressed(stream, T=arrays["T"], H_lm=arrays["H_lm"], taper=arrays["taper"],
                                metadata=json.dumps(changed))
        with pytest.raises(ContractError): load_ffi_result(tampered, allow_test=True)

    nonfinite = arrays["H_lm"].copy(); nonfinite[3] = np.nan
    bad = tmp_path / "nonfinite.npz"
    with bad.open("wb") as stream:
        np.savez_compressed(stream, T=arrays["T"], H_lm=nonfinite, taper=arrays["taper"],
                            metadata=arrays["metadata"])
    with pytest.raises(ContractError): load_ffi_result(bad, allow_test=True)

    old_metadata = dict(original_metadata); old_metadata.pop("transform_start")
    old = tmp_path / "old_conflated_metadata.npz"
    with old.open("wb") as stream:
        np.savez_compressed(stream, T=arrays["T"], H_lm=arrays["H_lm"], taper=arrays["taper"],
                            metadata=json.dumps(old_metadata))
    with pytest.raises(ContractError): load_ffi_result(old, allow_test=True)

    bad_taper = arrays["taper"].copy(); bad_taper[10] = 0.123
    bad = tmp_path / "bad_taper.npz"
    changed = dict(original_metadata)
    changed["output_hash"] = hashlib.sha256(arrays["T"].tobytes() + arrays["H_lm"].tobytes()
                                                + bad_taper.tobytes()).hexdigest()
    with bad.open("wb") as stream:
        np.savez_compressed(stream, T=arrays["T"], H_lm=arrays["H_lm"], taper=bad_taper,
                            metadata=json.dumps(changed))
    with pytest.raises(ContractError): load_ffi_result(bad, allow_test=True)

    with pytest.raises(ContractError):
        write_ffi_result(tmp_path / "false_complete.npz", result, status="complete",
                         source_mode_metadata=source_metadata)
    wrong_mode = dict(source_metadata); wrong_mode["m"] = 1
    with pytest.raises(ContractError):
        write_ffi_result(tmp_path / "wrong_mode.npz", result, status="test",
                         source_mode_metadata=wrong_mode)

    changed = dict(original_metadata)
    changed["source_mode_metadata"] = {**source_metadata, "phi_dot_at_analysis_start": "invalid"}
    tampered = tmp_path / "tampered_source_numeric.npz"
    with tampered.open("wb") as stream:
        np.savez_compressed(stream, T=arrays["T"], H_lm=arrays["H_lm"], taper=arrays["taper"],
                            metadata=json.dumps(changed))
    with pytest.raises(ContractError):
        load_ffi_result(tampered, allow_test=True)
