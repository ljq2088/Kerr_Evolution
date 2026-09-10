from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
import numpy as np
from kerr_waveform_tools import quintic_source_turn_on

from .errors import ContractError
from .geometry import horizons
from .trajectory import Jet3, StageSample

ArrayDerivative = Callable[[np.ndarray], np.ndarray]


@dataclass(frozen=True)
class SourceBlocks:
    block1: np.ndarray
    block2: np.ndarray
    block3: np.ndarray
    block4: np.ndarray
    horizon_total_limit: np.ndarray | None = None

    @property
    def total(self) -> np.ndarray:
        result = self.block1 + self.block2 + self.block3 + self.block4
        if self.horizon_total_limit is not None:
            mask = np.isfinite(self.horizon_total_limit)
            result = np.where(mask, self.horizon_total_limit, result)
        return result

    @classmethod
    def zeros(cls, shape: tuple[int, ...]) -> "SourceBlocks":
        items = [np.zeros(shape, dtype=complex) for _ in range(4)]
        return cls(*items)


@dataclass(frozen=True)
class SourceGeometry:
    R: np.ndarray
    y: np.ndarray
    zero_R: np.ndarray
    horizon: np.ndarray
    radial_horizon: float
    c: np.ndarray
    ell0: np.ndarray
    ellm1: np.ndarray
    jT: np.ndarray
    jR: np.ndarray
    u1_factor: np.ndarray
    w1_factor: np.ndarray
    h1_factor: np.ndarray
    u2_factor: np.ndarray
    w2_factor: np.ndarray
    h2_factor: np.ndarray
    u3_factor: np.ndarray
    h3_factor: np.ndarray
    u4_factor: np.ndarray
    w4_factor: np.ndarray
    prefactor: np.ndarray


def _regrouped_horizon_total(
    geometry: SourceGeometry,
    h1: np.ndarray,
    h2: np.ndarray,
    h3: np.ndarray,
    h4: np.ndarray,
) -> np.ndarray:
    direct_total = geometry.prefactor * (h1 + h2 + h3 + h4)
    result = np.full(geometry.R.shape, np.nan, dtype=complex)
    result[geometry.horizon] = direct_total[geometry.horizon]
    if not np.all(np.isfinite(result[geometry.horizon])):
        raise FloatingPointError("non-finite direct A5 total horizon limit")
    return result


def _source_geometry(R, y, M: float, a: float, L: float, m: int) -> SourceGeometry:
    R, y = np.broadcast_arrays(np.asarray(R, dtype=float), np.asarray(y, dtype=float))
    radial_horizon = L**2 / horizons(M, a)[0]
    zero_R = R == 0
    safe_R = np.where(zero_R, 1.0, R)
    rho = safe_R / (L**2 + 1j * a * safe_R * y)
    rhobar = safe_R / (L**2 - 1j * a * safe_R * y)
    c = np.sqrt(1 - y**2)
    rho8_over_R4 = np.where(zero_R, 0.0, safe_R**4 / (L**2 + 1j * a * safe_R * y) ** 8)
    h2_factor = -rho8_over_R4 * rhobar / np.sqrt(2)
    return SourceGeometry(
        R=R,
        y=y,
        zero_R=zero_R,
        horizon=np.isclose(R, radial_horizon, rtol=0.0, atol=8 * np.finfo(float).eps),
        radial_horizon=radial_horizon,
        c=c,
        ell0=m / c,
        ellm1=(m + y) / c,
        jT=-(2 + 4 * M * R / L**2),
        jR=-R**2 / L**2,
        u1_factor=rho**-2 * rhobar**-1,
        w1_factor=rho**-4,
        h1_factor=-rho**8 * rhobar,
        u2_factor=rho**-2 * rhobar**-2 * R**2,
        w2_factor=rho**-4 * rhobar**2,
        h2_factor=h2_factor,
        u3_factor=rho**-2 * rhobar,
        h3_factor=-rho8_over_R4 * rhobar / 2,
        u4_factor=rho**-2 * rhobar**-2,
        w4_factor=rho**-4 * rhobar**2 * R**2,
        prefactor=-16 * np.pi * R * (L**4 + a**2 * R**2 * y**2),
    )


def evaluate_regrouped_horizon_total_limit(
    R, y, sample: StageSample, M: float, a: float, L: float, m: int,
    sigma_R: float, sigma_y: float, dR: ArrayDerivative, dy: ArrayDerivative,
    epsilon_g: float,
    *,
    _geometry: SourceGeometry | None = None,
    _dressed: tuple[tuple[np.ndarray, np.ndarray, np.ndarray], ...] | None = None,
) -> np.ndarray:
    """Direct A5 total regular limit at R=R_H, independent of generic block outputs."""
    geometry = _source_geometry(R, y, M, a, L, m) if _geometry is None else _geometry
    R, y, horizon = geometry.R, geometry.y, geometry.horizon
    if not np.any(horizon):
        raise ContractError("horizon-limit evaluator requires an R=R_H row")
    fnn, fmn, fmm = (_all_dressed_jets(geometry, sample, m, sigma_R, sigma_y, epsilon_g)
                      if _dressed is None else _dressed)

    u10, u11, u12 = (geometry.u1_factor * value for value in fnn)
    v10 = geometry.c * dy(u10) - 1j * a * geometry.c * u11 + geometry.ell0 * u10
    v11 = geometry.c * dy(u11) - 1j * a * geometry.c * u12 + geometry.ell0 * u11
    w10, w11 = geometry.w1_factor * v10, geometry.w1_factor * v11
    h1 = geometry.h1_factor * (
        geometry.c * dy(w10) - 1j * a * geometry.c * w11 + geometry.ellm1 * w10)

    u20, u21, u22 = (geometry.u2_factor * value for value in fmn)
    v20 = geometry.jT * u21 + geometry.jR * dR(u20)
    v21 = geometry.jT * u22 + geometry.jR * dR(u21)
    w20, w21 = geometry.w2_factor * v20, geometry.w2_factor * v21
    h2 = geometry.h2_factor * (
        geometry.c * dy(w20) - 1j * a * geometry.c * w21 + geometry.ellm1 * w20)

    u30, u31, u32 = (geometry.u3_factor * value for value in fmm)
    v30 = geometry.jT * u31 + geometry.jR * dR(u30)
    v31 = geometry.jT * u32 + geometry.jR * dR(u31)
    w30, w31 = geometry.w1_factor * v30, geometry.w1_factor * v31
    h3 = geometry.h3_factor * (geometry.jT * w31 + geometry.jR * dR(w30))

    u40, u41, u42 = (geometry.u4_factor * value for value in fmn)
    v40 = geometry.c * dy(u40) - 1j * a * geometry.c * u41 + geometry.ellm1 * u40
    v41 = geometry.c * dy(u41) - 1j * a * geometry.c * u42 + geometry.ellm1 * u41
    w40, w41 = geometry.w4_factor * v40, geometry.w4_factor * v41
    h4 = geometry.h2_factor * (geometry.jT * w41 + geometry.jR * dR(w40))

    return _regrouped_horizon_total(geometry, h1, h2, h3, h4)


def quintic_window(T: float, tau_on: float = 20.0) -> Jet3:
    jet = quintic_source_turn_on(T, tau_on)
    return Jet3(jet.value, jet.first, jet.second)


def window_amplitude(amplitude: Jet3, window: Jet3) -> Jet3:
    return Jet3(window.value * amplitude.value,
                window.first * amplitude.value + window.value * amplitude.first,
                window.second * amplitude.value + 2 * window.first * amplitude.first + window.value * amplitude.second)


def _phase_jet(amplitude: Jet3, phase: Jet3, m: int) -> Jet3:
    factor = np.exp(-1j * m * phase.value)
    return Jet3(
        factor * amplitude.value,
        factor * (amplitude.first - 1j * m * phase.first * amplitude.value),
        factor * (amplitude.second - 2j * m * phase.first * amplitude.first - 1j * m * phase.second * amplitude.value - m**2 * phase.first**2 * amplitude.value),
    )


def gaussian_dressed_jets(R, y, sample: StageSample, amplitude: Jet3, m: int, sigma_R: float, sigma_y: float, epsilon_g: float = 0.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    R, y = np.broadcast_arrays(np.asarray(R, dtype=float), np.asarray(y, dtype=float))
    gaussian, gamma1, gamma2 = _gaussian_kinematics(R, y, sample, sigma_R, sigma_y, epsilon_g)
    return _dress_amplitude(amplitude, sample, m, gaussian, gamma1, gamma2)


def _gaussian_kinematics(R, y, sample: StageSample, sigma_R: float, sigma_y: float, epsilon_g: float):
    if sigma_R <= 0 or sigma_y <= 0:
        raise ContractError("Gaussian widths must be positive")
    x = R - sample.R
    exponent = np.exp(-0.5 * (x**2 / sigma_R**2 + y**2 / sigma_y**2))
    if epsilon_g:
        if not 0 < epsilon_g < 1:
            raise ContractError("epsilon_g must lie in (0,1)")
        exponent = np.where(exponent >= epsilon_g, exponent, 0.0)
    gaussian = exponent / (2 * np.pi * sigma_R * sigma_y)
    gamma1 = x * sample.R_jet.first / sigma_R**2
    gamma2 = x**2 * sample.R_jet.first**2 / sigma_R**4 + (x * sample.R_jet.second - sample.R_jet.first**2) / sigma_R**2
    return gaussian, gamma1, gamma2


def _dress_amplitude(amplitude: Jet3, sample: StageSample, m: int, gaussian, gamma1, gamma2):
    q = _phase_jet(amplitude, sample.Phi_jet, m)
    return q.value * gaussian, (q.first + q.value * gamma1) * gaussian, (q.second + 2 * q.first * gamma1 + q.value * gamma2) * gaussian


def _all_dressed_jets(geometry: SourceGeometry, sample: StageSample, m: int, sigma_R: float, sigma_y: float, epsilon_g: float):
    kinematics = _gaussian_kinematics(
        geometry.R, geometry.y, sample, sigma_R, sigma_y, epsilon_g)
    return tuple(
        _dress_amplitude(amplitude, sample, m, *kinematics)
        for amplitude in (sample.amplitude_nn, sample.amplitude_mn, sample.amplitude_mm)
    )


@dataclass(frozen=True)
class _LocalField:
    data: np.ndarray
    r0: int
    y0: int
    full_shape: tuple[int, int]

    @property
    def r1(self) -> int:
        return self.r0 + self.data.shape[0]

    @property
    def y1(self) -> int:
        return self.y0 + self.data.shape[1]

    def multiply(self, factor) -> "_LocalField":
        if np.ndim(factor) == 0:
            data = self.data * factor
        else:
            data = self.data * factor[self.r0:self.r1, self.y0:self.y1]
        return _LocalField(data, self.r0, self.y0, self.full_shape)

    def to_full(self) -> np.ndarray:
        result = np.zeros(self.full_shape, dtype=complex)
        result[self.r0:self.r1, self.y0:self.y1] = self.data
        return result


def _add_local(*fields: _LocalField) -> _LocalField:
    r0, r1 = min(item.r0 for item in fields), max(item.r1 for item in fields)
    y0, y1 = min(item.y0 for item in fields), max(item.y1 for item in fields)
    result = np.zeros((r1 - r0, y1 - y0), dtype=complex)
    for item in fields:
        result[item.r0 - r0:item.r1 - r0, item.y0 - y0:item.y1 - y0] += item.data
    return _LocalField(result, r0, y0, fields[0].full_shape)


class _LocalSourceDerivatives:
    def __init__(self, grid) -> None:
        self.grid = grid
        self._radial = {}
        self._angular = {}

    @staticmethod
    def _restriction(matrix, start: int, stop: int):
        rows = np.flatnonzero(np.asarray(matrix[:, start:stop].getnnz(axis=1)).ravel())
        if not len(rows) or not np.array_equal(rows, np.arange(rows[0], rows[-1] + 1)):
            raise ContractError("local source stencil support must be nonempty and contiguous")
        first, final = int(rows[0]), int(rows[-1]) + 1
        return first, matrix[first:final]

    def dr(self, field: _LocalField) -> _LocalField:
        key = (field.r0, field.r1)
        first, matrix = self._radial.setdefault(
            key, self._restriction(self.grid.D1R, *key))
        padded = np.zeros((field.full_shape[0], field.data.shape[1]), dtype=complex)
        padded[field.r0:field.r1] = field.data
        return _LocalField(matrix @ padded, first, field.y0, field.full_shape)

    def dy(self, field: _LocalField) -> _LocalField:
        key = (field.y0, field.y1)
        first, matrix = self._angular.setdefault(
            key, self._restriction(self.grid.D1y_source, *key))
        padded = np.zeros((field.full_shape[1], field.data.shape[0]), dtype=complex)
        padded[field.y0:field.y1] = field.data.T
        return _LocalField((matrix @ padded).T, field.r0, first, field.full_shape)


def _local_dressed_jets(
    geometry: SourceGeometry,
    sample: StageSample,
    m: int,
    sigma_R: float,
    sigma_y: float,
    epsilon_g: float,
):
    if not 0 < epsilon_g < 1:
        raise ContractError("local source evaluation requires epsilon_g in (0,1)")
    radial = np.exp(-0.5 * ((geometry.R[:, 0] - sample.R) / sigma_R) ** 2) >= epsilon_g
    angular = np.exp(-0.5 * (geometry.y[0] / sigma_y) ** 2) >= epsilon_g
    radial_indices, angular_indices = np.flatnonzero(radial), np.flatnonzero(angular)
    if not len(radial_indices) or not len(angular_indices):
        return None
    r0, r1 = int(radial_indices[0]), int(radial_indices[-1]) + 1
    y0, y1 = int(angular_indices[0]), int(angular_indices[-1]) + 1
    R = geometry.R[r0:r1, y0:y1]
    y = geometry.y[r0:r1, y0:y1]
    gaussian, gamma1, gamma2 = _gaussian_kinematics(
        R, y, sample, sigma_R, sigma_y, epsilon_g)
    shape = geometry.R.shape
    return tuple(
        tuple(_LocalField(value, r0, y0, shape) for value in
              _dress_amplitude(amplitude, sample, m, gaussian, gamma1, gamma2))
        for amplitude in (sample.amplitude_nn, sample.amplitude_mn, sample.amplitude_mm)
    )


def _evaluate_local_source_blocks(
    geometry: SourceGeometry,
    derivatives: _LocalSourceDerivatives,
    sample: StageSample,
    a: float,
    m: int,
    sigma_R: float,
    sigma_y: float,
    epsilon_g: float,
    zero_blocks: SourceBlocks,
) -> SourceBlocks:
    cutoff = np.sqrt(2 * np.log(1 / epsilon_g))
    if np.any(geometry.zero_R) and sample.R <= cutoff * sigma_R:
        raise ContractError("retained Gaussian support reaches SCRI+; analytic zero source row is unavailable")
    dressed = _local_dressed_jets(geometry, sample, m, sigma_R, sigma_y, epsilon_g)
    if dressed is None:
        return zero_blocks
    f_nn, f_mn, f_mm = dressed

    u10, u11, u12 = (item.multiply(geometry.u1_factor) for item in f_nn)
    v10 = _add_local(
        derivatives.dy(u10).multiply(geometry.c),
        u11.multiply(-1j * a * geometry.c),
        u10.multiply(geometry.ell0),
    )
    v11 = _add_local(
        derivatives.dy(u11).multiply(geometry.c),
        u12.multiply(-1j * a * geometry.c),
        u11.multiply(geometry.ell0),
    )
    w10, w11 = v10.multiply(geometry.w1_factor), v11.multiply(geometry.w1_factor)
    h1 = _add_local(
        derivatives.dy(w10).multiply(geometry.c),
        w11.multiply(-1j * a * geometry.c),
        w10.multiply(geometry.ellm1),
    ).multiply(geometry.h1_factor)

    u20, u21, u22 = (item.multiply(geometry.u2_factor) for item in f_mn)
    v20 = _add_local(
        u21.multiply(geometry.jT),
        derivatives.dr(u20).multiply(geometry.jR),
    )
    v21 = _add_local(
        u22.multiply(geometry.jT),
        derivatives.dr(u21).multiply(geometry.jR),
    )
    w20, w21 = v20.multiply(geometry.w2_factor), v21.multiply(geometry.w2_factor)
    h2 = _add_local(
        derivatives.dy(w20).multiply(geometry.c),
        w21.multiply(-1j * a * geometry.c),
        w20.multiply(geometry.ellm1),
    ).multiply(geometry.h2_factor)

    u30, u31, u32 = (item.multiply(geometry.u3_factor) for item in f_mm)
    v30 = _add_local(
        u31.multiply(geometry.jT),
        derivatives.dr(u30).multiply(geometry.jR),
    )
    v31 = _add_local(
        u32.multiply(geometry.jT),
        derivatives.dr(u31).multiply(geometry.jR),
    )
    w30, w31 = v30.multiply(geometry.w1_factor), v31.multiply(geometry.w1_factor)
    h3 = _add_local(
        w31.multiply(geometry.jT),
        derivatives.dr(w30).multiply(geometry.jR),
    ).multiply(geometry.h3_factor)

    u40, u41, u42 = (item.multiply(geometry.u4_factor) for item in f_mn)
    v40 = _add_local(
        derivatives.dy(u40).multiply(geometry.c),
        u41.multiply(-1j * a * geometry.c),
        u40.multiply(geometry.ellm1),
    )
    v41 = _add_local(
        derivatives.dy(u41).multiply(geometry.c),
        u42.multiply(-1j * a * geometry.c),
        u41.multiply(geometry.ellm1),
    )
    w40, w41 = v40.multiply(geometry.w4_factor), v41.multiply(geometry.w4_factor)
    h4 = _add_local(
        w41.multiply(geometry.jT),
        derivatives.dr(w40).multiply(geometry.jR),
    ).multiply(geometry.h2_factor)

    blocks = tuple(item.multiply(geometry.prefactor).to_full() for item in (h1, h2, h3, h4))
    for block in blocks:
        block[geometry.zero_R] = 0.0
    regrouped = _add_local(h1, h2, h3, h4).multiply(geometry.prefactor).to_full()
    horizon_limit = np.full(geometry.R.shape, np.nan, dtype=complex)
    horizon_limit[geometry.horizon] = regrouped[geometry.horizon]
    if not np.all(np.isfinite(horizon_limit[geometry.horizon])):
        raise FloatingPointError("non-finite local A5 total horizon limit")
    return SourceBlocks(*blocks, horizon_limit)


def evaluate_source_blocks(
    R,
    y,
    sample: StageSample,
    M: float,
    a: float,
    L: float,
    m: int,
    sigma_R: float,
    sigma_y: float,
    dR: ArrayDerivative,
    dy: ArrayDerivative,
    *,
    epsilon_g: float = 0.0,
    allow_boundary_limits: bool = False,
    _geometry: SourceGeometry | None = None,
) -> SourceBlocks:
    geometry = _source_geometry(R, y, M, a, L, m) if _geometry is None else _geometry
    R, y = geometry.R, geometry.y
    if (not allow_boundary_limits and (np.any(R <= 0) or np.any(R >= geometry.radial_horizon))) or np.any(R < 0) or np.any(R > geometry.radial_horizon) or np.any(np.abs(y) >= 1):
        raise ContractError("A5 four-block evaluator is restricted to 0 < R < R_H and |y| < 1")
    if np.any(geometry.zero_R) and sample.R <= np.sqrt(2 * np.log(1 / max(epsilon_g, np.finfo(float).tiny))) * sigma_R:
        raise ContractError("retained Gaussian support reaches SCRI+; analytic zero source row is unavailable")
    f_nn, f_mn, f_mm = _all_dressed_jets(geometry, sample, m, sigma_R, sigma_y, epsilon_g)

    u10, u11, u12 = (geometry.u1_factor * item for item in f_nn)
    v10 = geometry.c * dy(u10) - 1j * a * geometry.c * u11 + geometry.ell0 * u10
    v11 = geometry.c * dy(u11) - 1j * a * geometry.c * u12 + geometry.ell0 * u11
    w10, w11 = geometry.w1_factor * v10, geometry.w1_factor * v11
    h1 = geometry.h1_factor * (
        geometry.c * dy(w10) - 1j * a * geometry.c * w11 + geometry.ellm1 * w10)

    u20, u21, u22 = (geometry.u2_factor * item for item in f_mn)
    v20 = geometry.jT * u21 + geometry.jR * dR(u20)
    v21 = geometry.jT * u22 + geometry.jR * dR(u21)
    w20, w21 = geometry.w2_factor * v20, geometry.w2_factor * v21
    h2 = geometry.h2_factor * (
        geometry.c * dy(w20) - 1j * a * geometry.c * w21 + geometry.ellm1 * w20)

    u30, u31, u32 = (geometry.u3_factor * item for item in f_mm)
    v30 = geometry.jT * u31 + geometry.jR * dR(u30)
    v31 = geometry.jT * u32 + geometry.jR * dR(u31)
    w30, w31 = geometry.w1_factor * v30, geometry.w1_factor * v31
    h3 = geometry.h3_factor * (geometry.jT * w31 + geometry.jR * dR(w30))

    u40, u41, u42 = (geometry.u4_factor * item for item in f_mn)
    v40 = geometry.c * dy(u40) - 1j * a * geometry.c * u41 + geometry.ellm1 * u40
    v41 = geometry.c * dy(u41) - 1j * a * geometry.c * u42 + geometry.ellm1 * u41
    w40, w41 = geometry.w4_factor * v40, geometry.w4_factor * v41
    h4 = geometry.h2_factor * (geometry.jT * w41 + geometry.jR * dR(w40))
    blocks = [np.asarray(geometry.prefactor * h, dtype=complex) for h in (h1, h2, h3, h4)]
    for block in blocks:
        block[geometry.zero_R] = 0.0
    horizon_limit = (_regrouped_horizon_total(geometry, h1, h2, h3, h4)
                     if allow_boundary_limits and np.any(geometry.horizon) else None)
    return SourceBlocks(*blocks, horizon_limit)


class GridSourceEvaluator:
    def __init__(self, grid, trajectory, m: int, sigma_R: float, sigma_y: float, epsilon_g: float = 1e-14, tau_on: float = 20.0):
        if not 0 < tau_on <= 40:
            raise ContractError("startup time must satisfy 0 < tau_on <= 40M")
        self.grid, self.trajectory, self.m = grid, trajectory, m
        self.sigma_R, self.sigma_y, self.epsilon_g, self.tau_on = sigma_R, sigma_y, epsilon_g, tau_on
        self._geometry = _source_geometry(grid.RR, grid.yy, 1.0, trajectory.a, 1.0, m)
        self._local_derivatives = _LocalSourceDerivatives(grid)
        cutoff_cells_R = int(np.ceil(np.sqrt(2 * np.log(1 / epsilon_g)) * sigma_R / grid.dR))
        cutoff_cells_y = int(np.ceil(np.sqrt(2 * np.log(1 / epsilon_g)) * sigma_y / grid.dy))
        local_rows = min(len(grid.R), 2 * cutoff_cells_R + 13)
        local_columns = min(len(grid.y), 2 * cutoff_cells_y + 13)
        self.uses_local_support = local_rows * local_columns < 0.45 * grid.RR.size
        self._zero_blocks = SourceBlocks.zeros(grid.RR.shape)
        self._zero_total = np.zeros(grid.RR.shape, dtype=complex)
        for array in (*self._zero_blocks.__dict__.values(), self._zero_total):
            if isinstance(array, np.ndarray):
                array.setflags(write=False)

    def blocks(self, T: float) -> SourceBlocks:
        if T >= self.trajectory.events.source_off:
            return self._zero_blocks
        sample = self.trajectory.sample(T)
        if sample.R >= self.trajectory.R_horizon:
            return self._zero_blocks
        window = quintic_window(T, self.tau_on)
        sample = StageSample(sample.T, sample.r, sample.R, sample.Phi, sample.R_jet, sample.Phi_jet,
                             sample.u_T, sample.u_R, sample.u_Phi,
                             window_amplitude(sample.amplitude_nn, window),
                             window_amplitude(sample.amplitude_mn, window),
                             window_amplitude(sample.amplitude_mm, window))
        if self.uses_local_support:
            return _evaluate_local_source_blocks(
                self._geometry, self._local_derivatives, sample, self.trajectory.a, self.m,
                self.sigma_R, self.sigma_y, self.epsilon_g, self._zero_blocks)
        return evaluate_source_blocks(
            self.grid.RR, self.grid.yy, sample, 1.0, self.trajectory.a, 1.0, self.m,
            self.sigma_R, self.sigma_y, self.grid.dr, self.grid.dy_source,
            epsilon_g=self.epsilon_g, allow_boundary_limits=True, _geometry=self._geometry)

    def __call__(self, T: float) -> np.ndarray:
        if T >= self.trajectory.events.source_off:
            return self._zero_total
        return self.blocks(T).total
