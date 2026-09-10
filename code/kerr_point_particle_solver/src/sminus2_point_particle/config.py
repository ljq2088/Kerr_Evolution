from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json

from kerr_waveform_tools.contracts import ContractError
from kerr_waveform_tools.transition import (
    TRANSITION_LMAX,
    TRANSITION_MASS_RATIO,
    TRANSITION_PRESCRIPTION,
    TRANSITION_T_OT,
    TRANSITION_X_START,
    TransitionSpec,
)


INTERNAL_MASS = 1.0
COMPACTIFICATION_SCALE = 1.0
T0 = 0.0
PHI0 = 0.0


@dataclass(frozen=True)
class EvolutionConfig:
    mass_scale: float
    chi: float
    m: int
    n_r: int
    n_y: int
    post_light_ring_time_over_m: float = 120.0
    sigma_r_over_dr: float = 4.0
    sigma_y_over_dy: float = 4.0
    transition_prescription: str = TRANSITION_PRESCRIPTION
    flux_model_version: str = "unassigned"
    tau_on_over_m: float = 20.0
    gaussian_tail_tolerance: float = 1e-14

    def __post_init__(self) -> None:
        if self.mass_scale <= 0 or not 0 <= self.chi < 1:
            raise ContractError("mass_scale must be positive and 0 <= chi < 1")
        if self.n_r < 3 or self.n_y < 3:
            raise ContractError("grids need >= 3 nodes")
        if self.m not in (2, 4):
            raise ContractError("current factor-aware implementation supports only target modes m=2 or m=4")
        if self.post_light_ring_time_over_m <= 0:
            raise ContractError("post-light-ring duration must be explicitly positive")
        if self.sigma_r_over_dr <= 0 or self.sigma_y_over_dy <= 0:
            raise ContractError("Gaussian width ratios must be positive")
        if not 0 < self.tau_on_over_m <= 40:
            raise ContractError("startup tau_on must satisfy 0 < tau_on <= 40M")
        if not 0 < self.gaussian_tail_tolerance < 1:
            raise ContractError("Gaussian tail tolerance must lie in (0,1)")
        if self.transition_prescription != TRANSITION_PRESCRIPTION:
            raise ContractError(f"C0 fixes transition_prescription={TRANSITION_PRESCRIPTION!r}")

    @property
    def M(self) -> float:
        return INTERNAL_MASS

    @property
    def a(self) -> float:
        return self.chi

    @property
    def L(self) -> float:
        return COMPACTIFICATION_SCALE

    @property
    def q_mass(self) -> float:
        return TRANSITION_MASS_RATIO

    def to_json(self) -> str:
        record = asdict(self)
        record["q_mass"] = TRANSITION_MASS_RATIO
        return json.dumps(record, sort_keys=True, separators=(",", ":"))

    @classmethod
    def from_json(cls, payload: str) -> "EvolutionConfig":
        record = json.loads(payload)
        q_mass = record.pop("q_mass", TRANSITION_MASS_RATIO)
        if q_mass != TRANSITION_MASS_RATIO:
            raise ContractError(f"C0 fixes transition q_mass={TRANSITION_MASS_RATIO!r}")
        return cls(**record)

    def initial_data_key(self, schema_version: str, generator_version: str) -> str:
        return TransitionSpec(
            self.chi,
            self.flux_model_version,
            M=INTERNAL_MASS,
            L=COMPACTIFICATION_SCALE,
            transition_prescription=self.transition_prescription,
        ).initial_data_key(schema_version, generator_version)


@dataclass(frozen=True)
class OutputRequest:
    ell_out: tuple[int, ...] = (2,)
    save_field: bool = False
    save_scri_field: bool = False
    save_psi4_lm: bool = True
    save_strain_lm: bool = False
    output_dt_over_m: float = 0.1
    field_dt_over_m: float = 5.0
    checkpoint_dt_over_m: float = 10.0
    chunk_size: int = 1024

    def __post_init__(self) -> None:
        if not self.ell_out or any(ell < 2 for ell in self.ell_out):
            raise ContractError("spin -2 modes require ell >= 2")
        if min(self.output_dt_over_m, self.field_dt_over_m, self.checkpoint_dt_over_m, self.chunk_size) <= 0:
            raise ContractError("cadences and chunk_size must be positive")
        for name, cadence in (("field", self.field_dt_over_m), ("checkpoint", self.checkpoint_dt_over_m)):
            if not abs(cadence / self.output_dt_over_m - round(cadence / self.output_dt_over_m)) < 1e-12:
                raise ContractError(f"{name} cadence must be an integer multiple of output cadence")

    def to_json(self) -> str:
        record = asdict(self)
        record["ell_out"] = list(self.ell_out)
        return json.dumps(record, sort_keys=True, separators=(",", ":"))

    @classmethod
    def from_json(cls, payload: str) -> "OutputRequest":
        record = json.loads(payload)
        record["ell_out"] = tuple(record["ell_out"])
        return cls(**record)
