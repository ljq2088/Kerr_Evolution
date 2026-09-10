from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import os
import hashlib
from pathlib import Path
from typing import Protocol
import numpy as np

from .contracts import ContractError
from .background import horizons, isco_quantities, radial_state


TRANSITION_PRESCRIPTION = "ori-thorne-v1"
TRANSITION_MASS_RATIO = 1.0e-5
TRANSITION_LMAX = 6
TRANSITION_T_OT = 3.412
TRANSITION_X_START = -0.5381356808642278
TRANSITION_SCHEMA_VERSION = "1"
TRANSITION_GENERATOR_VERSION = "ori-thorne-python-v1"


@dataclass(frozen=True)
class TransitionSpec:
    chi: float
    flux_model_version: str
    M: float = 1.0
    L: float = 1.0
    transition_prescription: str = TRANSITION_PRESCRIPTION

    def __post_init__(self) -> None:
        if (not 0 <= self.chi < 1 or self.M <= 0 or self.L <= 0
                or self.flux_model_version in ("", "unassigned")
                or self.transition_prescription != TRANSITION_PRESCRIPTION):
            raise ContractError("invalid Kerr transition specification")

    @property
    def a(self) -> float:
        return self.chi * self.M

    def initial_data_key(self, schema_version: str, generator_version: str) -> str:
        if not schema_version or not generator_version:
            raise ContractError("schema and generator versions are required for an initial-data key")
        record = {
            "chi": self.chi,
            "M": self.M,
            "L": self.L,
            "background_hash": background_hash(self.chi, self.M, self.L),
            "q_mass": TRANSITION_MASS_RATIO,
            "lmax": TRANSITION_LMAX,
            "t_ot": TRANSITION_T_OT,
            "x_start": TRANSITION_X_START,
            "transition_prescription": self.transition_prescription,
            "flux_model_version": self.flux_model_version,
            "schema_version": schema_version,
            "generator_version": generator_version,
        }
        body = json.dumps(record, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(body.encode("ascii")).hexdigest()


@dataclass(frozen=True)
class TransitionManifest:
    prescription: str
    flux_model_version: str
    delta_r: float
    delta_energy: float
    delta_lz: float
    q_mass: float = TRANSITION_MASS_RATIO
    lmax: int = TRANSITION_LMAX
    t_ot: float = TRANSITION_T_OT
    x_start: float = TRANSITION_X_START
    flux_total_over_mu2: float | None = None
    relative_lmax_increment: float | None = None
    schema_version: str = TRANSITION_SCHEMA_VERSION
    generator_version: str = TRANSITION_GENERATOR_VERSION
    chi: float | None = None
    M: float = 1.0
    L: float = 1.0
    background_hash: str = ""
    r_isco: float | None = None
    energy_isco: float | None = None
    angular_momentum_isco: float | None = None
    omega_isco: float | None = None
    edot_correction: float | None = None
    alpha_ot: float | None = None
    beta_ot: float | None = None
    kappa_ot: float | None = None
    tau0: float | None = None
    radial_scale: float | None = None
    flux_provenance: str = ""


@dataclass(frozen=True)
class PlungeInitialData:
    r0: float
    energy: float
    angular_momentum: float
    manifest: TransitionManifest

    def as_record(self) -> dict:
        return asdict(self)


class InitialDataCache(Protocol):
    def load(self, key: str) -> PlungeInitialData | None: ...
    def store(self, key: str, value: PlungeInitialData) -> None: ...


def background_hash(chi: float, M: float, L: float) -> str:
    payload = json.dumps({"chi": chi, "M": M, "L": L}, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def validate_initial_data(config, value: PlungeInitialData) -> PlungeInitialData:
    manifest = value.manifest
    if manifest.chi != config.chi or manifest.M != config.M or manifest.L != config.L:
        raise ContractError("initial data background does not match EvolutionConfig")
    if manifest.flux_model_version != config.flux_model_version:
        raise ContractError("initial data flux version does not match EvolutionConfig")
    expected = initial_data_from_manifest(config.M, config.a, config.L, manifest)
    if value != expected:
        raise ContractError("stored plunge values do not match manifest-derived values")
    return value


class MemoryInitialDataCache:
    def __init__(self, config) -> None:
        self.config = config
        self._records: dict[str, PlungeInitialData] = {}

    def load(self, key: str) -> PlungeInitialData | None:
        value = self._records.get(key)
        if value is None:
            return None
        self._validate_pair(key, value)
        return value

    def store(self, key: str, value: PlungeInitialData) -> None:
        self._validate_pair(key, value)
        self._records[key] = value

    def _validate_pair(self, key: str, value: PlungeInitialData) -> None:
        validate_initial_data(self.config, value)
        expected = self.config.initial_data_key(value.manifest.schema_version, value.manifest.generator_version)
        if key != expected:
            raise ContractError("cache key does not match config and manifest versions")


class JsonInitialDataCache:
    """Versioned, atomic one-record-per-key cache."""

    def __init__(self, directory: Path, config) -> None:
        self.directory = Path(directory)
        self.config = config

    def _path(self, key: str) -> Path:
        if len(key) != 64 or any(c not in "0123456789abcdef" for c in key):
            raise ContractError("cache key must be a lowercase SHA-256 digest")
        return self.directory / f"{key}.json"

    def load(self, key: str) -> PlungeInitialData | None:
        path = self._path(key)
        if not path.exists():
            return None
        record = json.loads(path.read_text(encoding="ascii"))
        record["manifest"] = TransitionManifest(**record["manifest"])
        value = PlungeInitialData(**record)
        self._validate_pair(key, value)
        return value

    def store(self, key: str, value: PlungeInitialData) -> None:
        self._validate_pair(key, value)
        self.directory.mkdir(parents=True, exist_ok=True)
        path = self._path(key)
        if path.exists():
            if self.load(key) != value:
                raise ContractError("refusing to overwrite a different initial-data record")
            return
        temporary = path.with_suffix(".tmp")
        temporary.write_text(json.dumps(value.as_record(), sort_keys=True, indent=2) + "\n", encoding="ascii")
        os.replace(temporary, path)

    def _validate_pair(self, key: str, value: PlungeInitialData) -> None:
        validate_initial_data(self.config, value)
        expected = self.config.initial_data_key(value.manifest.schema_version, value.manifest.generator_version)
        if key != expected:
            raise ContractError("cache key does not match config and manifest versions")


@dataclass(frozen=True)
class FluxRecord:
    chi: float
    lmax: int
    flux_total_over_mu2: float
    relative_lmax_increment: float
    flux_version: str
    provenance: str = "unspecified"

    def __post_init__(self) -> None:
        if not 0 <= self.chi < 1 or self.lmax != TRANSITION_LMAX:
            raise ContractError(f"flux input requires 0 <= chi < 1 and lmax={TRANSITION_LMAX}")
        if (self.flux_total_over_mu2 <= 0 or self.relative_lmax_increment < 0 or not self.flux_version
                or self.provenance in ("", "unspecified")):
            raise ContractError("flux values must be positive/nonnegative and explicitly versioned")


def generate_ori_thorne_initial_data(flux: FluxRecord, M: float = 1.0, L: float = 1.0) -> PlungeInitialData:
    if M != 1.0 or L != 1.0:
        raise ContractError("the target-consistent generator currently supports internal M=L=1 only")
    r, energy, angular_momentum, omega = isco_quantities(flux.chi, M)
    a = flux.chi * M
    edot = (5 / 32) * flux.flux_total_over_mu2 * (M * omega) ** (-10 / 3)
    alpha = 3 / r**6 * (r**2 + 2 * (a**2 * (energy**2 - 1) - angular_momentum**2) * r + 10 * (angular_momentum - a * energy) ** 2)
    beta = 2 / r**4 * ((angular_momentum - a**2 * energy * omega) * r - 3 * (angular_momentum - a * energy) * (1 - a * omega))
    kappa = (32 / 5) * omega ** (7 / 3) * (1 + a / r**1.5) / (1 - 3 * M / r + 2 * a * np.sqrt(M) / r**1.5) ** 0.5 * edot
    if min(alpha, beta, kappa) <= 0:
        raise ContractError("Ori--Thorne scales require positive alpha, beta, and kappa")
    tau0 = (alpha * beta * kappa) ** (-1 / 5)
    radial_scale = (beta * kappa) ** (2 / 5) * alpha ** (-3 / 5)
    delta_lz = kappa * tau0 * TRANSITION_T_OT * TRANSITION_MASS_RATIO ** (4 / 5)
    manifest = TransitionManifest(
        TRANSITION_PRESCRIPTION, flux.flux_version,
        abs(TRANSITION_X_START) * TRANSITION_MASS_RATIO ** (2 / 5) * radial_scale,
        omega * delta_lz, delta_lz,
        flux_total_over_mu2=flux.flux_total_over_mu2,
        relative_lmax_increment=flux.relative_lmax_increment,
        chi=flux.chi, M=M, L=L, background_hash=background_hash(flux.chi, M, L),
        r_isco=r, energy_isco=energy, angular_momentum_isco=angular_momentum, omega_isco=omega,
        edot_correction=edot, alpha_ot=alpha, beta_ot=beta, kappa_ot=kappa,
        tau0=tau0, radial_scale=radial_scale, flux_provenance=flux.provenance,
    )
    return initial_data_from_manifest(M, a, L, manifest)


def initial_data_from_manifest(M: float, a: float, L: float, manifest: TransitionManifest) -> PlungeInitialData:
    chi = a / M
    if manifest.chi != chi or manifest.M != M or manifest.L != L or manifest.background_hash != background_hash(chi, M, L):
        raise ContractError("transition manifest background fields/hash do not match requested Kerr background")
    if manifest.prescription != TRANSITION_PRESCRIPTION:
        raise ContractError(f"C0 accepts only the fixed {TRANSITION_PRESCRIPTION!r} manifest")
    if manifest.q_mass != TRANSITION_MASS_RATIO:
        raise ContractError(f"C0 fixes transition q_mass={TRANSITION_MASS_RATIO!r}")
    if (manifest.lmax, manifest.t_ot, manifest.x_start) != (TRANSITION_LMAX, TRANSITION_T_OT, TRANSITION_X_START):
        raise ContractError("transition manifest does not match the fixed target-consistent constants")
    if manifest.flux_model_version in ("", "unassigned"):
        raise ContractError("transition manifest requires a versioned flux model")
    if (manifest.schema_version != TRANSITION_SCHEMA_VERSION
            or manifest.generator_version != TRANSITION_GENERATOR_VERSION):
        raise ContractError("unsupported transition schema or generator version")
    required = (manifest.r_isco, manifest.energy_isco, manifest.angular_momentum_isco, manifest.omega_isco,
                manifest.edot_correction, manifest.alpha_ot, manifest.beta_ot, manifest.kappa_ot,
                manifest.tau0, manifest.radial_scale, manifest.flux_total_over_mu2,
                manifest.relative_lmax_increment)
    if any(item is None or not np.isfinite(item) for item in required) or not manifest.flux_provenance:
        raise ContractError("transition manifest is missing complete ISCO/OT/flux provenance")
    r_isco, e_isco, lz_isco, omega_isco = isco_quantities(chi, M)
    if not np.allclose((manifest.r_isco, manifest.energy_isco, manifest.angular_momentum_isco, manifest.omega_isco),
                       (r_isco, e_isco, lz_isco, omega_isco), rtol=2e-14, atol=2e-14):
        raise ContractError("manifest ISCO quantities do not match its background")
    edot = (5 / 32) * manifest.flux_total_over_mu2 * (M * omega_isco) ** (-10 / 3)
    alpha = 3 / r_isco**6 * (r_isco**2 + 2 * (a**2 * (e_isco**2 - 1) - lz_isco**2) * r_isco
                              + 10 * (lz_isco - a * e_isco) ** 2)
    beta = 2 / r_isco**4 * ((lz_isco - a**2 * e_isco * omega_isco) * r_isco
                             - 3 * (lz_isco - a * e_isco) * (1 - a * omega_isco))
    kappa = (32 / 5) * omega_isco ** (7 / 3) * (1 + a / r_isco**1.5) / (
        1 - 3 * M / r_isco + 2 * a * np.sqrt(M) / r_isco**1.5) ** 0.5 * edot
    tau0 = (alpha * beta * kappa) ** (-1 / 5)
    radial_scale = (beta * kappa) ** (2 / 5) * alpha ** (-3 / 5)
    delta_lz = kappa * tau0 * TRANSITION_T_OT * TRANSITION_MASS_RATIO ** (4 / 5)
    derived = (edot, alpha, beta, kappa, tau0, radial_scale,
               abs(TRANSITION_X_START) * TRANSITION_MASS_RATIO ** (2 / 5) * radial_scale,
               omega_isco * delta_lz, delta_lz)
    recorded = (manifest.edot_correction, manifest.alpha_ot, manifest.beta_ot, manifest.kappa_ot,
                manifest.tau0, manifest.radial_scale, manifest.delta_r, manifest.delta_energy, manifest.delta_lz)
    if not np.allclose(recorded, derived, rtol=3e-14, atol=3e-14):
        raise ContractError("manifest OT scales or deltas do not recompute from its versioned flux")
    value = PlungeInitialData(r_isco - manifest.delta_r, e_isco - manifest.delta_energy, lz_isco - manifest.delta_lz, manifest)
    rp, _ = horizons(M, a)
    if not rp < value.r0 < r_isco:
        raise ContractError("transition manifest must yield r_+ < r0 < r_ISCO")
    state = radial_state(value.r0, M, a, L, value.energy, value.angular_momentum)
    omega_h = a / (2 * M * rp)
    if state.potential < 0 or value.energy - omega_h * value.angular_momentum <= 0:
        raise ContractError("transition manifest violates A6 plunge inequalities")
    return value
