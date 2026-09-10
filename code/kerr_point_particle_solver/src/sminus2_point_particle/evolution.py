from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.optimize import minimize_scalar

from .errors import ContractError
from .field import FieldCoefficients, angular_operator, coefficients, pi_from_p


@dataclass(frozen=True)
class ContinuumRHS:
    dP: np.ndarray
    dpsi: np.ndarray


def continuum_rhs(psi, P, Q, Q_R, angular_psi, source, coeff: FieldCoefficients) -> ContinuumRHS:
    dpsi = pi_from_p(psi, P, Q, coeff)
    dP = source - coeff.C * Q_R + angular_psi - coeff.E_R * Q - coeff.F * psi
    return ContinuumRHS(np.asarray(dP, dtype=complex), np.asarray(dpsi, dtype=complex))


@dataclass(frozen=True)
class EvolutionState:
    P: np.ndarray
    psi: np.ndarray


@dataclass(frozen=True)
class TimeStepChoice:
    dt: float
    steps_per_output: int
    dt_wave: float
    dt_motion: float
    dt_phase: float
    dt_amplitude: float
    controlling_limit: str
    q_reference: float
    q_floor: float
    gamma_source: float
    t_stable: float
    extrema_interval: tuple[float, float]
    dt_ramp: float


class SpatialRHS:
    def __init__(self, grid, chi: float, m: int, source=None):
        self.grid, self.m, self.source = grid, m, source
        self.coeff = coefficients(grid.RR, grid.yy, 1.0, chi, 1.0, m)

    def __call__(self, T: float, state: EvolutionState) -> EvolutionState:
        Q = self.grid.dr(state.psi)
        Q_R = self.grid.drr(state.psi)
        angular = angular_operator(state.psi, self.grid.dy_field(state.psi), self.grid.dyy_field(state.psi), self.grid.yy, self.m)
        source = np.zeros_like(state.psi) if self.source is None else self.source(T)
        result = continuum_rhs(state.psi, state.P, Q, Q_R, angular, source, self.coeff)
        if not np.all(np.isfinite(result.dP)) or not np.all(np.isfinite(result.dpsi)):
            raise FloatingPointError(f"non-finite RHS at T={T:.16g}")
        return EvolutionState(result.dP, result.dpsi)


def zero_state(shape: tuple[int, int]) -> EvolutionState:
    return EvolutionState(np.zeros(shape, dtype=complex), np.zeros(shape, dtype=complex))


def rk4_step(rhs, T: float, state: EvolutionState, dt: float, *, stage_observer=None) -> EvolutionState:
    def stage(time, candidate):
        if stage_observer is not None:
            stage_observer(time)
        return rhs(time, candidate)
    def shifted(base, rate, scale):
        return EvolutionState(base.P + scale * rate.P, base.psi + scale * rate.psi)
    k1 = stage(T, state)
    k2 = stage(T + dt / 2, shifted(state, k1, dt / 2))
    k3 = stage(T + dt / 2, shifted(state, k2, dt / 2))
    k4 = stage(T + dt, shifted(state, k3, dt))
    return EvolutionState(state.P + dt * (k1.P + 2 * k2.P + 2 * k3.P + k4.P) / 6,
                          state.psi + dt * (k1.psi + 2 * k2.psi + 2 * k3.psi + k4.psi) / 6)


def choose_timestep(grid, chi: float, m: int, trajectory, output_dt: float = 0.1, *, tau_on: float = 20.0,
                    default_t_stable: float = 100.0, sigma_R: float | None = None) -> TimeStepChoice:
    if not 0 < tau_on <= 40:
        raise ContractError("timestep startup time must satisfy 0 < tau_on <= 40M")
    coeff = coefficients(grid.RR, grid.yy, 1.0, chi, 1.0, m)
    discriminant = np.maximum(np.real(coeff.B**2 / 4 - coeff.A * coeff.C), 0.0)
    speeds = [np.abs((coeff.B / 2 + sign * np.sqrt(discriminant)) / coeff.A) for sign in (1, -1)]
    vr = max(float(np.max(speed)) for speed in speeds)
    vy = float(np.max(np.sqrt((1 - grid.yy**2) / coeff.A)))
    dt_wave = 0.20 / np.sqrt((vr / grid.dR) ** 2 + (vy / grid.dy) ** 2)
    horizon_time = trajectory.events.horizon_crossing
    if horizon_time <= tau_on:
        raise ContractError("horizon crossing must occur after source startup completes")
    t_stable = max(tau_on, horizon_time - 20.0) if horizon_time < 120.0 else default_t_stable
    if not tau_on <= t_stable < horizon_time:
        raise ContractError("no nonempty fully-on interval is available for global source timestep extrema")
    nodes = np.asarray(trajectory.T)
    nodes = nodes[(nodes > t_stable) & (nodes < horizon_time)]
    nodes = np.unique(np.concatenate(([t_stable], nodes, [np.nextafter(horizon_time, -np.inf)])))

    def conservative_max(function) -> float:
        intervals = list(zip(nodes[:-1], nodes[1:]))
        scores = []
        maximum = 0.0
        for left, right in intervals:
            values = [function(left), function((left + right) / 2), function(right)]
            score = max(values); scores.append(score); maximum = max(maximum, score)
        for index in np.argsort(scores)[-min(16, len(scores)):]:
            left, right = intervals[int(index)]
            result = minimize_scalar(lambda T: -function(float(T)), bounds=(left, right), method="bounded",
                                     options={"xatol": 1e-12})
            maximum = max(maximum, -float(result.fun))
        return 1.01 * maximum

    max_rdot = conservative_max(lambda T: abs(trajectory.sample(float(T)).R_jet.first))
    max_phidot = conservative_max(lambda T: abs(trajectory.sample(float(T)).Phi_jet.first))
    sigma_R = 4 * grid.dR if sigma_R is None else sigma_R
    dt_motion = 0.50 * sigma_R / max_rdot
    dt_phase = np.inf if m == 0 else 0.10 / (abs(m) * max_phidot)
    def q_pairs(T):
        sample = trajectory.sample(float(T))
        phase = np.exp(-1j * m * sample.Phi)
        pairs = []
        for raw in (sample.amplitude_nn, sample.amplitude_mn, sample.amplitude_mm):
            pairs.append((phase * raw.value,
                          phase * (raw.first - 1j * m * sample.Phi_jet.first * raw.value)))
        return pairs

    q_reference = conservative_max(lambda T: max(abs(pair[0]) for pair in q_pairs(T)))
    q_floor = 1e-12 * q_reference
    gamma = 0.0 if q_reference == 0 else conservative_max(
        lambda T: max(abs(q1) / max(abs(q0), q_floor) for q0, q1 in q_pairs(T)))
    dt_amplitude = np.inf if gamma == 0 else 0.15 / gamma
    dt_ramp = tau_on / 40.0
    limits = {"wave": dt_wave, "motion": dt_motion, "phase": dt_phase,
              "amplitude": dt_amplitude, "ramp": dt_ramp}
    controlling = min(limits, key=limits.get)
    steps = int(np.ceil(output_dt / limits[controlling]))
    return TimeStepChoice(output_dt / steps, steps, dt_wave, dt_motion, dt_phase, dt_amplitude,
                          controlling, q_reference, q_floor, gamma, t_stable,
                          (t_stable, horizon_time), dt_ramp)


def evolve_accepted_steps(rhs, initial: EvolutionState, start: float, stop: float, dt: float, *, observer=None) -> EvolutionState:
    span = (stop - start) / dt
    steps = int(round(span))
    if dt <= 0 or not np.isclose(span, steps, rtol=0, atol=1e-11):
        raise ContractError("evolution interval must contain an integer number of RK4 steps")
    state, T = initial, float(start)
    for _ in range(steps):
        state = rk4_step(rhs, T, state, dt)
        T += dt
        if observer is not None:
            observer(T, state)
    return state
