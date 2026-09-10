from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol
import numpy as np
from scipy.integrate import solve_ivp

from .contracts import ContractError
from .background import RadialState, compactified_radius, horizons, radial_state
from .transition import PlungeInitialData, validate_initial_data


@dataclass(frozen=True)
class Jet3:
    value: complex
    first: complex
    second: complex


@dataclass(frozen=True)
class StageSample:
    T: float
    r: float
    R: float
    Phi: float
    R_jet: Jet3
    Phi_jet: Jet3
    u_T: Jet3
    u_R: Jet3
    u_Phi: Jet3
    amplitude_nn: Jet3
    amplitude_mn: Jet3
    amplitude_mm: Jet3


class DenseTrajectory(Protocol):
    def sample(self, T: float) -> StageSample: ...


def _product_jet(A0: Jet3, X: Jet3, Y: Jet3) -> Jet3:
    return Jet3(
        A0.value * X.value * Y.value,
        A0.first * X.value * Y.value + A0.value * (X.first * Y.value + X.value * Y.first),
        A0.second * X.value * Y.value
        + 2 * A0.first * (X.first * Y.value + X.value * Y.first)
        + A0.value * (X.second * Y.value + 2 * X.first * Y.first + X.value * Y.second),
    )


def stage_sample(T: float, r: float, Phi: float, M: float, a: float, L: float, energy: float, angular_momentum: float, particle_mass: float = 1.0, *, allow_horizon_interior: bool = False) -> StageSample:
    if particle_mass < 0:
        raise ContractError("particle_mass must be nonnegative")
    state: RadialState = radial_state(r, M, a, L, energy, angular_momentum, allow_horizon_interior=allow_horizon_interior)
    rdot = -state.potential_sqrt / (r**2 * state.u_T)
    vprime = -(1 / L**2) * (
        (2 * r * state.u_R + r**2 * state.u_R_r) / state.u_T
        - r**2 * state.u_R * state.u_T_r / state.u_T**2
    )
    rddot = rdot * vprime
    R = float(compactified_radius(r, L))
    Rdot = state.u_R / state.u_T
    Rddot = rdot * (state.u_R_r * state.u_T - state.u_R * state.u_T_r) / state.u_T**2
    Phidot = state.u_Phi / state.u_T
    Phiddot = rdot * (state.u_Phi_r * state.u_T - state.u_Phi * state.u_T_r) / state.u_T**2
    def velocity_jet(value, first, second):
        return Jet3(value, first * rdot, second * rdot**2 + first * rddot)
    uT = velocity_jet(state.u_T, state.u_T_r, state.u_T_rr)
    uR = velocity_jet(state.u_R, state.u_R_r, state.u_R_rr)
    uP = velocity_jet(state.u_Phi, state.u_Phi_r, state.u_Phi_rr)
    nhat = -state.K / (2 * L**4)
    nhat1 = -state.K_r * rdot / (2 * L**4)
    nhat2 = -(state.K_rr * rdot**2 + state.K_r * rddot) / (2 * L**4)
    B = angular_momentum - a * energy
    mbar = -1j * B / (np.sqrt(2) * r)
    mbar1 = 1j * B * rdot / (np.sqrt(2) * r**2)
    mbar2 = -2j * B * rdot**2 / (np.sqrt(2) * r**3) + 1j * B * rddot / (np.sqrt(2) * r**2)
    A0 = particle_mass * L**2 / (2 * np.pi * r**4 * state.u_T)
    lam = 4 * rdot / r + uT.first / state.u_T
    lamdot = 4 * (rddot / r - rdot**2 / r**2) + uT.second / state.u_T - (uT.first / state.u_T) ** 2
    weight = Jet3(A0, -lam * A0, (lam**2 - lamdot) * A0)
    nh = Jet3(nhat, nhat1, nhat2)
    mb = Jet3(mbar, mbar1, mbar2)
    return StageSample(
        T, r, R, Phi, Jet3(R, Rdot, Rddot), Jet3(Phi, Phidot, Phiddot), uT, uR, uP,
        _product_jet(weight, nh, nh), _product_jet(weight, mb, nh), _product_jet(weight, mb, mb),
    )


def trajectory_rates(R: float, Phi: float, M: float, a: float, L: float, energy: float, angular_momentum: float) -> np.ndarray:
    if R <= 0:
        raise ContractError("trajectory rates require the finite-radius domain R > 0")
    r = L**2 / R
    state = radial_state(r, M, a, L, energy, angular_momentum, allow_horizon_interior=True)
    return np.array([state.u_R / state.u_T, state.u_Phi / state.u_T], dtype=float)


def prograde_light_ring_radius(chi: float, M: float = 1.0) -> float:
    return 2 * M * (1 + np.cos((2 / 3) * np.arccos(-chi)))


@dataclass(frozen=True)
class TrajectoryEvents:
    light_ring: float
    horizon_crossing: float
    source_off: float
    terminal: float


@dataclass(frozen=True)
class DOP853Trajectory:
    T: np.ndarray
    R: np.ndarray
    Phi: np.ndarray
    events: TrajectoryEvents
    R_horizon: float
    R_terminal: float
    dense_solution: object
    M: float
    a: float
    L: float
    initial_data: PlungeInitialData
    particle_mass: float
    light_ring_scri_reference: float

    def position(self, T: float) -> tuple[float, float]:
        if T < self.T[0] or T > self.events.terminal:
            raise ContractError("stage time lies outside the retained dense trajectory")
        return tuple(map(float, self.dense_solution(float(T))))

    def sample(self, T: float) -> StageSample:
        R, Phi = self.position(T)
        return stage_sample(T, self.L**2 / R, Phi, self.M, self.a, self.L,
                            self.initial_data.energy, self.initial_data.angular_momentum,
                            self.particle_mass, allow_horizon_interior=True)

    def evolution_end_time(self, post_light_ring: float = 120.0) -> float:
        return self.light_ring_scri_reference + post_light_ring


def tortoise_radius(r: float, M: float, a: float) -> float:
    rp, rm = horizons(M, a)
    return (r + 2 * M * rp / (rp - rm) * np.log(abs((r - rp) / M))
            - 2 * M * rm / (rp - rm) * np.log(abs((r - rm) / M)))


def integrate_trajectory(
    initial_data: PlungeInitialData,
    config,
    sigma_R: float,
    *,
    particle_mass: float = 1.0,
    rtol: float = 2e-11,
    atol: float = 2e-13,
    max_step: float = 0.25,
    max_time: float = 1000.0,
) -> DOP853Trajectory:
    validate_initial_data(config, initial_data)
    chi = config.chi
    M = L = 1.0
    a = chi
    r_plus = horizons(M, a)[0]
    R_horizon = L**2 / r_plus
    R_lr = L**2 / prograde_light_ring_radius(chi, M)
    R_terminal = R_horizon + np.sqrt(2 * np.log(1e14)) * sigma_R

    def rhs(_T, state):
        return trajectory_rates(float(state[0]), float(state[1]), M, a, L, initial_data.energy, initial_data.angular_momentum)

    def event_at(target, terminal=False):
        def event(_T, state):
            return state[0] - target
        event.direction = 1.0
        event.terminal = terminal
        return event

    solution = solve_ivp(rhs, (0.0, max_time), (L**2 / initial_data.r0, 0.0), method="DOP853",
                         rtol=rtol, atol=atol, max_step=max_step, dense_output=True,
                         events=(event_at(R_lr), event_at(R_horizon), event_at(R_terminal, True)))
    if not solution.success or solution.sol is None or any(len(events) != 1 for events in solution.t_events):
        raise ContractError(f"DOP853 trajectory failed to record each required event once: {solution.message}")
    event_times = [float(events[0]) for events in solution.t_events]
    if not event_times[0] < event_times[1] < event_times[2]:
        raise ContractError("trajectory events are not ordered light-ring < horizon < terminal")
    if not np.all(np.diff(solution.y[0]) > 0) or not np.all(np.isfinite(solution.y)):
        raise ContractError("trajectory is not finite and monotonically increasing in R")
    events = TrajectoryEvents(event_times[0], event_times[1], event_times[1], event_times[2])
    r_lr = prograde_light_ring_radius(chi, M)
    rstar_lr = tortoise_radius(r_lr, M, a)
    height_lr = rstar_lr - 2 * r_lr - 4 * M * np.log(r_lr / M)
    scri_reference = event_times[0] - height_lr - rstar_lr
    return DOP853Trajectory(np.asarray(solution.t), np.asarray(solution.y[0]), np.unwrap(solution.y[1]), events,
                            R_horizon, R_terminal, solution.sol, M, a, L, initial_data, particle_mass, scri_reference)
