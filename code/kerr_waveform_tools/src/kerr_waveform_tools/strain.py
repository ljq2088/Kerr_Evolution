from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import sys
import numpy as np

from .artifacts import ArtifactStatus, nonnegative_count
from .contracts import ContractError


FFI_SCHEMA = "c2-H-v2"
SUPPORTED_FFI_SCHEMAS = {"c2-H-v1", FFI_SCHEMA}
FFT_CONVENTION = "numpy fft: forward exp(-2pi i k n/N), inverse 1/N"
SAMPLING_STATUSES = {"owner-waived", "executed"}
CUTOFF_DIAGNOSTIC_STATUSES = {"not_run", "executed"}


@dataclass(frozen=True)
class FFIMetadata:
    ell: int
    m: int
    dt: float
    omega_i: float
    omega0: float
    taper_width: float
    L2_bridge: float
    a1_status: str
    input_hash: str
    transform_start: float
    transform_stop: float
    analysis_start: float
    analysis_stop: float
    fft_convention: str
    cutoff_factor: float
    sampling_status: str
    sampling_epsilon: float | None
    sampling_threshold: float | None
    input_schema: str
    source_input_hash: str
    source_events: dict
    source_normalization: str
    source_error_metadata: dict
    source_environment: dict
    source_background: dict
    amplitude_scaling: str
    recommended_start: float
    recommended_stop: float
    cutoff_diagnostic_status: str
    transform_interval: tuple[float, float]
    analysis_interval: tuple[float, float]
    recommended_untapered_interval: tuple[float, float]
    implementation_hash: str
    implementation_environment: dict
    output_object: str = "H_lm=r*h_lm"


@dataclass(frozen=True)
class FFIResult:
    T: np.ndarray
    H: np.ndarray
    taper: np.ndarray
    metadata: FFIMetadata


def _validated_uniform_series(T, values, *, name: str) -> tuple[np.ndarray, np.ndarray, float]:
    T = np.asarray(T, dtype=float)
    values = np.asarray(values, dtype=complex)
    if T.ndim != 1 or values.ndim != 1 or T.shape != values.shape or len(T) < 4:
        raise ContractError(f"{name} requires equal one-dimensional arrays with at least four samples")
    if not np.all(np.isfinite(T)) or not np.all(np.isfinite(values.real)) or not np.all(np.isfinite(values.imag)):
        raise ContractError(f"{name} arrays must be finite")
    differences = np.diff(T)
    if np.any(differences <= 0):
        raise ContractError(f"{name} time must be strictly increasing")
    if not np.allclose(differences, differences[0], rtol=1e-12, atol=1e-14):
        raise ContractError(f"{name} time must be uniformly sampled")
    return T, values, float(differences[0])


def _implementation_provenance() -> tuple[str, dict]:
    digest = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    environment = {"python": sys.version.split()[0], "numpy": np.__version__}
    return digest, environment


def quintic_taper(T: np.ndarray, start: float, stop: float, width: float = 10.0) -> np.ndarray:
    T = np.asarray(T, dtype=float)
    if (T.ndim != 1 or not np.all(np.isfinite(T)) or not np.isfinite(start) or not np.isfinite(stop)
            or not np.isfinite(width) or width <= 0 or stop - start < 2 * width
            or len(T) == 0 or start < T[0] or stop > T[-1] or start >= stop):
        raise ContractError("FFI interval must accommodate two positive-width tapers")
    window = np.zeros_like(T)
    middle = (T >= start + width) & (T <= stop - width)
    window[middle] = 1
    left = (T >= start) & (T < start + width)
    x = (T[left] - start) / width
    window[left] = 10 * x**3 - 15 * x**4 + 6 * x**5
    right = (T > stop - width) & (T <= stop)
    x = (stop - T[right]) / width
    window[right] = 10 * x**3 - 15 * x**4 + 6 * x**5
    return window


def fixed_frequency_spectrum(frequencies, psi4_spectrum, *, omega0: float,
                             L: float = 1.0) -> np.ndarray:
    frequencies = np.asarray(frequencies, dtype=float)
    psi4_spectrum = np.asarray(psi4_spectrum, dtype=complex)
    if (
        frequencies.ndim != 1
        or psi4_spectrum.shape != frequencies.shape
        or not np.all(np.isfinite(frequencies))
        or not np.all(np.isfinite(psi4_spectrum))
        or not np.isfinite(omega0)
        or omega0 <= 0
        or not np.isfinite(L)
        or L <= 0
    ):
        raise ContractError("invalid direct fixed-frequency spectrum input")
    return np.asarray(
        -(L**2) * psi4_spectrum / np.maximum(np.abs(frequencies), omega0) ** 2,
        dtype=np.complex128,
    )


def fixed_frequency_integrate(T, psi4_lm, *, ell: int, m: int, phi_dot_at_start: float,
                              transform_start: float | None = None, analysis_start: float = 100.0,
                              stop: float | None = None, L: float = 1.0,
                              cutoff_factor: float = 0.75, taper_width: float = 10.0,
                              a1_status: str = "conditional-open", input_schema: str = "array-api",
                              source_input_hash: str = "unavailable", source_events: dict | None = None,
                              source_normalization: str = "unspecified", source_error_metadata: dict | None = None,
                              source_environment: dict | None = None,
                              source_background: dict | None = None,
                              amplitude_scaling: str = "unspecified",
                              sampling_status: str = "owner-waived", sampling_epsilon: float | None = None,
                              sampling_threshold: float | None = None,
                              cutoff_diagnostic_status: str = "not_run") -> FFIResult:
    T, psi4_lm, dt = _validated_uniform_series(T, psi4_lm, name="FFI")
    transform_start = float(T[0] if transform_start is None else transform_start)
    analysis_start = float(analysis_start)
    stop = float(T[-1] if stop is None else stop)
    if (isinstance(ell, (bool, np.bool_)) or isinstance(m, (bool, np.bool_))
            or not isinstance(ell, (int, np.integer)) or not isinstance(m, (int, np.integer)) or ell < 2
            or abs(m) > ell or not np.isfinite(phi_dot_at_start) or not np.isfinite(L) or L <= 0
            or not np.isfinite(cutoff_factor) or cutoff_factor <= 0):
        raise ContractError("invalid FFI mode, trajectory frequency, L, or cutoff factor")
    if a1_status != "conditional-open":
        raise ContractError("C2 requires the unresolved A1 bridge to remain conditional-open")
    if sampling_status not in SAMPLING_STATUSES:
        raise ContractError("sampling status must be owner-waived or executed")
    if cutoff_diagnostic_status not in CUTOFF_DIAGNOSTIC_STATUSES:
        raise ContractError("cutoff diagnostic status must be not_run or executed")
    if sampling_status == "owner-waived" and (sampling_epsilon is not None or sampling_threshold is not None):
        raise ContractError("owner-waived sampling cannot carry executed diagnostic values")
    if sampling_status == "executed":
        values = (sampling_epsilon, sampling_threshold)
        if (any(value is None or isinstance(value, (bool, np.bool_)) or not np.isfinite(value) for value in values)
                or sampling_epsilon < 0 or sampling_threshold <= 0):
            raise ContractError("executed sampling requires finite nonnegative epsilon and positive threshold")
    if (not np.isclose(transform_start, T[0], rtol=0.0, atol=1e-12)
            or not transform_start + taper_width <= analysis_start <= stop - taper_width):
        raise ContractError("transform must start at the first real sample and contain the analysis interval")
    omega_i = abs(m * phi_dot_at_start)
    omega0 = cutoff_factor * omega_i
    if omega0 <= 0:
        raise ContractError("FFI cutoff requires nonzero |m*Phi_dot(T_a)|")
    window = quintic_taper(T, transform_start, stop, taper_width)
    omega = 2 * np.pi * np.fft.fftfreq(len(T), d=dt)
    denominator = np.maximum(np.abs(omega), omega0) ** 2
    H = np.fft.ifft(-(L**2) * np.fft.fft(window * psi4_lm) / denominator)
    digest = hashlib.sha256(T.tobytes() + psi4_lm.tobytes()).hexdigest()
    implementation_hash, implementation_environment = _implementation_provenance()
    metadata = FFIMetadata(int(ell), int(m), dt, omega_i, omega0, taper_width, L**2, a1_status, digest,
                           transform_start, stop, analysis_start, stop, FFT_CONVENTION,
                           cutoff_factor, sampling_status, sampling_epsilon, sampling_threshold,
                           input_schema, source_input_hash, source_events or {}, source_normalization,
                           source_error_metadata or {}, source_environment or {}, source_background or {},
                           amplitude_scaling, analysis_start,
                           stop - taper_width, cutoff_diagnostic_status, (transform_start, stop),
                           (analysis_start, stop), (analysis_start, stop - taper_width),
                           implementation_hash, implementation_environment)
    return FFIResult(T, H, window, metadata)


def write_ffi_result(path: Path, result: FFIResult, *, status: str, source_mode_metadata: dict) -> None:
    try:
        artifact_status = ArtifactStatus(status)
    except (ValueError, TypeError) as error:
        raise ContractError("invalid C2 output status") from error
    if artifact_status not in (ArtifactStatus.TEST, ArtifactStatus.COMPLETE):
        raise ContractError("C2 output must be test or complete")
    required_source = {"schema", "status", "ell", "m", "events", "normalization", "error_metadata",
                       "input_hash", "environment", "run_kind", "T1", "analysis_start",
                       "phi_dot_at_analysis_start", "background", "amplitude_scaling"}
    if not isinstance(source_mode_metadata, dict):
        raise ContractError("C2 source metadata must be a mapping")
    if missing := required_source - source_mode_metadata.keys():
        raise ContractError(f"C2 source metadata missing keys: {sorted(missing)}")
    try:
        source_status = ArtifactStatus(source_mode_metadata["status"])
    except (ValueError, TypeError) as error:
        raise ContractError("C2 source metadata has unknown status") from error
    if artifact_status is ArtifactStatus.COMPLETE and (
            source_status is not ArtifactStatus.COMPLETE or source_mode_metadata["run_kind"] != "full"):
        raise ContractError("complete C2 output requires a complete full-run mode")
    if (result.metadata.ell != source_mode_metadata["ell"]
            or result.metadata.m != source_mode_metadata["m"]):
        raise ContractError("C2 result/source mode mismatch")
    for value, key in ((result.metadata.input_schema, "schema"),
                       (result.metadata.source_input_hash, "input_hash"),
                       (result.metadata.source_events, "events"),
                       (result.metadata.source_normalization, "normalization"),
                       (result.metadata.source_error_metadata, "error_metadata"),
                       (result.metadata.source_environment, "environment"),
                       (result.metadata.source_background, "background"),
                       (result.metadata.amplitude_scaling, "amplitude_scaling")):
        if value != source_mode_metadata[key]:
            raise ContractError(f"C2 result/source provenance mismatch for {key}")
    T, H, _ = _validated_uniform_series(result.T, result.H, name="C2 output")
    taper = np.asarray(result.taper, dtype=float)
    if taper.shape != T.shape or not np.all(np.isfinite(taper)) or np.any((taper < 0) | (taper > 1)):
        raise ContractError("C2 taper must be finite, in [0,1], and match T")
    if (not np.isclose(result.metadata.analysis_start, source_mode_metadata["analysis_start"],
                       rtol=0.0, atol=1e-12)
            or not np.isclose(result.metadata.analysis_stop, source_mode_metadata["T1"],
                              rtol=0.0, atol=1e-12)
            or not np.isclose(T[-1], source_mode_metadata["T1"], rtol=0.0, atol=1e-12)
            or not np.isclose(result.metadata.transform_start, T[0], rtol=0.0, atol=1e-12)
            or not np.isclose(
                result.metadata.omega_i,
                abs(result.metadata.m * source_mode_metadata["phi_dot_at_analysis_start"]),
                rtol=1e-12, atol=1e-14)):
        raise ContractError("C2 result/source time or frequency metadata mismatch")
    output_hash = hashlib.sha256(T.tobytes() + H.tobytes() + taper.tobytes()).hexdigest()
    metadata = {**result.metadata.__dict__, "schema": FFI_SCHEMA, "status": artifact_status.value,
                "sample_count": len(T), "output_hash": output_hash,
                "source_mode_metadata": source_mode_metadata}
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, T=T, H_lm=H, taper=taper, metadata=json.dumps(metadata, sort_keys=True))
    os.replace(temporary, path)


def load_ffi_result(path: Path, *, expected_ell: int | None = None, expected_m: int | None = None,
                    allow_test: bool = False) -> tuple[FFIResult, dict]:
    with np.load(path, allow_pickle=False) as data:
        T, H, taper = np.asarray(data["T"]), np.asarray(data["H_lm"]), np.asarray(data["taper"])
        metadata = json.loads(str(data["metadata"]))
    if metadata.get("schema") not in SUPPORTED_FFI_SCHEMAS:
        raise ContractError("C2 schema mismatch")
    try:
        status = ArtifactStatus(metadata.get("status"))
    except (ValueError, TypeError) as error:
        raise ContractError("unknown C2 status") from error
    if status is not ArtifactStatus.COMPLETE and not (allow_test and status is ArtifactStatus.TEST):
        raise ContractError("C2 loader refuses non-complete output")
    if T.dtype != np.float64 or H.dtype != np.complex128 or taper.dtype != np.float64:
        raise ContractError("C2 arrays must use float64/complex128/float64 dtypes")
    T, H, _ = _validated_uniform_series(T, H, name="C2 loader")
    if taper.shape != T.shape or not np.all(np.isfinite(taper)) or np.any((taper < 0) | (taper > 1)):
        raise ContractError("C2 taper is invalid")
    if nonnegative_count(metadata.get("sample_count"), "sample_count") != len(T):
        raise ContractError("C2 sample_count mismatch")
    expected_hash = hashlib.sha256(T.tobytes() + H.tobytes() + taper.tobytes()).hexdigest()
    if metadata.get("output_hash") != expected_hash:
        raise ContractError("C2 output hash mismatch")
    if ((expected_ell is not None and metadata.get("ell") != expected_ell)
            or (expected_m is not None and metadata.get("m") != expected_m)):
        raise ContractError("C2 mode mismatch")
    fields = {field.name for field in FFIMetadata.__dataclass_fields__.values()}
    if missing := fields - metadata.keys():
        raise ContractError(f"C2 metadata missing keys: {sorted(missing)}")
    source = metadata.get("source_mode_metadata")
    required_source = {"schema", "status", "ell", "m", "events", "normalization", "error_metadata",
                       "input_hash", "environment", "run_kind", "T1", "analysis_start",
                       "phi_dot_at_analysis_start", "background", "amplitude_scaling"}
    if (not isinstance(source, dict) or required_source - source.keys()
            or source.get("schema") != metadata["input_schema"]):
        raise ContractError("C2 source provenance schema is missing or inconsistent")
    try:
        source_status = ArtifactStatus(source["status"])
        source_T1 = float(source["T1"])
        source_start = float(source["analysis_start"])
        source_phi_dot = float(source["phi_dot_at_analysis_start"])
    except (ValueError, TypeError) as error:
        raise ContractError("C2 source provenance status or numeric fields are invalid") from error
    if not np.all(np.isfinite((source_T1, source_start, source_phi_dot))):
        raise ContractError("C2 source provenance numeric fields must be finite")
    for top, nested in (("source_input_hash", "input_hash"), ("source_events", "events"),
                        ("source_normalization", "normalization"),
                        ("source_error_metadata", "error_metadata"), ("source_environment", "environment"),
                        ("source_background", "background"), ("amplitude_scaling", "amplitude_scaling")):
        if metadata[top] != source.get(nested):
            raise ContractError(f"C2 source provenance mismatch for {top}")
    if metadata["ell"] != source.get("ell") or metadata["m"] != source.get("m"):
        raise ContractError("C2 source provenance mode mismatch")
    if status is ArtifactStatus.COMPLETE and (
            source_status is not ArtifactStatus.COMPLETE or source.get("run_kind") != "full"):
        raise ContractError("complete C2 output lacks complete full-run provenance")
    if status is ArtifactStatus.TEST and source_status not in (ArtifactStatus.TEST, ArtifactStatus.INCOMPLETE):
        raise ContractError("test C2 output has inconsistent source status")
    try:
        numeric_names = (
            "dt", "omega_i", "omega0", "taper_width", "L2_bridge", "transform_start", "transform_stop", "analysis_start",
            "analysis_stop", "cutoff_factor", "recommended_start", "recommended_stop")
        if any(isinstance(metadata[name], (bool, np.bool_)) or not isinstance(metadata[name], (int, float, np.number))
               for name in numeric_names):
            raise TypeError("numeric metadata must use numeric scalar types")
        numeric_values = {name: float(metadata[name]) for name in numeric_names}
    except (TypeError, ValueError) as error:
        raise ContractError("C2 numeric metadata are invalid") from error
    if not np.all(np.isfinite(list(numeric_values.values()))) or any(numeric_values[name] <= 0 for name in (
            "dt", "omega_i", "omega0", "taper_width", "L2_bridge", "cutoff_factor")):
        raise ContractError("C2 numeric metadata must be finite and positive")
    if (not np.isclose(metadata["transform_start"], T[0], rtol=0.0, atol=1e-12)
            or not np.isclose(metadata["transform_stop"], T[-1], rtol=0.0, atol=1e-12)
            or not metadata["transform_start"] + metadata["taper_width"] <= metadata["analysis_start"]
            or not metadata["analysis_start"] < metadata["analysis_stop"] <= T[-1]
            or not np.isclose(metadata["recommended_start"],
                              metadata["analysis_start"], rtol=0.0, atol=1e-12)
            or not np.isclose(metadata["recommended_stop"],
                              metadata["analysis_stop"] - metadata["taper_width"], rtol=0.0, atol=1e-12)
            or metadata["sampling_status"] not in SAMPLING_STATUSES
            or metadata["cutoff_diagnostic_status"] not in CUTOFF_DIAGNOSTIC_STATUSES
            or metadata["a1_status"] != "conditional-open"
            or not np.allclose(metadata["transform_interval"],
                               [metadata["transform_start"], metadata["transform_stop"]], rtol=0.0, atol=1e-12)
            or not np.allclose(metadata["analysis_interval"],
                               [metadata["analysis_start"], metadata["analysis_stop"]], rtol=0.0, atol=1e-12)
            or not np.allclose(metadata["recommended_untapered_interval"],
                               [metadata["recommended_start"], metadata["recommended_stop"]], rtol=0.0, atol=1e-12)
            or not np.isclose(metadata["dt"], T[1] - T[0], rtol=1e-12, atol=1e-14)
            or not np.isclose(metadata["omega0"], metadata["cutoff_factor"] * metadata["omega_i"],
                              rtol=1e-12, atol=1e-14)
            or metadata["fft_convention"] != FFT_CONVENTION
            or metadata["output_object"] != "H_lm=r*h_lm"
            or not np.isclose(metadata["analysis_start"], source_start, rtol=0.0, atol=1e-12)
            or not np.isclose(metadata["analysis_stop"], source_T1, rtol=0.0, atol=1e-12)
            or not np.isclose(
                metadata["omega_i"], abs(metadata["m"] * source_phi_dot),
                rtol=1e-12, atol=1e-14)):
        raise ContractError("C2 transform interval or status metadata is inconsistent")
    if metadata["sampling_status"] == "owner-waived":
        if metadata["sampling_epsilon"] is not None or metadata["sampling_threshold"] is not None:
            raise ContractError("owner-waived sampling carries executed diagnostic values")
    else:
        sampling_values = (metadata["sampling_epsilon"], metadata["sampling_threshold"])
        if (any(value is None or isinstance(value, bool) or not np.isfinite(value) for value in sampling_values)
                or metadata["sampling_epsilon"] < 0 or metadata["sampling_threshold"] <= 0):
            raise ContractError("executed sampling metadata are incomplete")
    expected_taper = quintic_taper(T, metadata["transform_start"], metadata["transform_stop"],
                                   metadata["taper_width"])
    if not np.array_equal(taper, expected_taper):
        raise ContractError("stored C2 taper does not match transform metadata")
    if not isinstance(metadata["input_hash"], str) or len(metadata["input_hash"]) != 64:
        raise ContractError("C2 input hash is invalid")
    if (not isinstance(metadata["implementation_hash"], str) or len(metadata["implementation_hash"]) != 64
            or not isinstance(metadata["implementation_environment"], dict)
            or not metadata["implementation_environment"].get("python")
            or not metadata["implementation_environment"].get("numpy")):
        raise ContractError("C2 implementation provenance is invalid")
    background = metadata["source_background"]
    required_background = {"mass_scale", "M_internal", "chi", "a_internal", "L_internal"}
    if (not isinstance(background, dict) or required_background - background.keys()
            or not np.all(np.isfinite([background[name] for name in required_background]))
            or background["M_internal"] <= 0 or background["L_internal"] <= 0
            or not 0 <= background["chi"] < 1
            or not np.isclose(background["a_internal"],
                              background["chi"] * background["M_internal"], rtol=1e-14, atol=1e-14)
            or not isinstance(metadata["amplitude_scaling"], str)
            or not metadata["amplitude_scaling"]):
        raise ContractError("C2 physical background or amplitude scaling provenance is invalid")
    ffi_payload = {name: metadata[name] for name in fields}
    ffi_payload["transform_interval"] = tuple(ffi_payload["transform_interval"])
    ffi_payload["analysis_interval"] = tuple(ffi_payload["analysis_interval"])
    ffi_payload["recommended_untapered_interval"] = tuple(ffi_payload["recommended_untapered_interval"])
    ffi_metadata = FFIMetadata(**ffi_payload)
    return FFIResult(T, H, taper, ffi_metadata), metadata


@dataclass(frozen=True)
class SamplingComparison:
    epsilon: float
    threshold: float
    passed: bool


@dataclass(frozen=True)
class CutoffVariation:
    factors: tuple[float, ...]
    relative_differences: tuple[float, ...]


def compare_cutoff_variation(T, psi4_lm, *, ell: int, m: int, phi_dot_at_start: float,
                             factors: tuple[float, ...] = (0.5, 0.75, 1.0), **ffi_kwargs) -> CutoffVariation:
    if 0.75 not in factors or any(not np.isfinite(value) or value <= 0 for value in factors):
        raise ContractError("cutoff diagnostic factors must be positive and include baseline 0.75")
    results = [fixed_frequency_integrate(T, psi4_lm, ell=ell, m=m, phi_dot_at_start=phi_dot_at_start,
                                         cutoff_factor=factor, cutoff_diagnostic_status="executed", **ffi_kwargs)
               for factor in factors]
    baseline_result = results[factors.index(0.75)]
    selected = ((baseline_result.T >= baseline_result.metadata.recommended_start)
                & (baseline_result.T <= baseline_result.metadata.recommended_stop))
    baseline = baseline_result.H[selected]
    norm = np.linalg.norm(baseline)
    if norm == 0:
        raise ContractError("cutoff diagnostic baseline has zero norm in the recommended interval")
    differences = tuple(float(np.linalg.norm(result.H[selected] - baseline) / norm) for result in results)
    return CutoffVariation(factors, differences)


def absolute_coarse_indices(T_fine, *, science_origin: float = 0.0, fine_dt: float = 0.05,
                            coarse_dt: float = 0.10, tolerance: float = 2e-11) -> np.ndarray:
    T_fine = np.asarray(T_fine, dtype=float)
    if T_fine.ndim != 1 or len(T_fine) < 2 or not np.allclose(np.diff(T_fine), fine_dt, rtol=0.0, atol=tolerance):
        raise ContractError("sampling comparison requires a uniform 0.05M fine grid")
    labels = np.rint((T_fine - science_origin) / coarse_dt).astype(int)
    targets = science_origin + labels * coarse_dt
    indices = np.flatnonzero(np.abs(T_fine - targets) <= tolerance)
    if len(indices) < 2 or not np.allclose(np.diff(T_fine[indices]), coarse_dt, rtol=0.0, atol=tolerance):
        raise ContractError("fine samples do not contain a uniform absolute 0.10M science subgrid")
    return indices


def compare_sampling(T_fine, psi4_fine, *, ell: int, m: int, phi_dot_at_start: float, epsilon_v2: float,
                     science_origin: float = 0.0, alignment_tolerance: float = 2e-11, **ffi_kwargs) -> SamplingComparison:
    if epsilon_v2 <= 0:
        raise ContractError("a positive V2 error estimate is required for the sampling gate")
    indices = absolute_coarse_indices(T_fine, science_origin=science_origin, tolerance=alignment_tolerance)
    fine = fixed_frequency_integrate(T_fine, psi4_fine, ell=ell, m=m, phi_dot_at_start=phi_dot_at_start, **ffi_kwargs)
    coarse = fixed_frequency_integrate(np.asarray(T_fine)[indices], np.asarray(psi4_fine)[indices], ell=ell, m=m,
                                       phi_dot_at_start=phi_dot_at_start, **ffi_kwargs)
    common_T = np.asarray(T_fine)[indices]
    selected = ((common_T >= fine.metadata.recommended_start)
                & (common_T <= fine.metadata.recommended_stop))
    reference = fine.H[indices][selected]
    if np.linalg.norm(reference) == 0:
        raise ContractError("sampling diagnostic reference has zero norm in the recommended interval")
    epsilon = float(np.linalg.norm(coarse.H[selected] - reference) / np.linalg.norm(reference))
    threshold = min(1e-3, 0.1 * epsilon_v2)
    return SamplingComparison(epsilon, threshold, epsilon <= threshold)
