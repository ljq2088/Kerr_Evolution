from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import numpy as np
from numpy.polynomial import Polynomial
from scipy import sparse

from .errors import ContractError
from .geometry import horizons
from .operators import finite_difference_weights


STENCIL_WIDTH = 7
GRID_PRESETS = {"smoke": (128, 33), "low": (256, 65), "medium": (512, 129), "high": (1024, 257)}


def pole_powers(m: int) -> tuple[int, int]:
    return abs(m - 2) // 2, abs(m + 2) // 2


def _ordinary_matrix(nodes: np.ndarray, derivative: int) -> sparse.csr_matrix:
    n = len(nodes)
    rows, columns_all, data = [], [], []
    half = STENCIL_WIDTH // 2
    for row in range(n):
        start = max(0, min(row - half, n - STENCIL_WIDTH))
        columns = np.arange(start, start + STENCIL_WIDTH)
        weights = finite_difference_weights(nodes[row], nodes[columns], derivative)
        rows.extend([row] * STENCIL_WIDTH); columns_all.extend(columns); data.extend(weights)
    return sparse.csr_matrix((data, (rows, columns_all)), shape=(n, n))


def _factor_weights(nodes: np.ndarray, x0: float, derivative: int, power: int, north: bool) -> tuple[np.ndarray, float]:
    factor = (Polynomial([1.0, 1.0]) if north else Polynomial([1.0, -1.0])) ** power
    scale = np.max(np.abs(nodes - x0))
    z = (nodes - x0) / scale
    offset = Polynomial([-x0 / scale, 1 / scale])
    matrix = np.vstack([factor(nodes) * z**k for k in range(STENCIL_WIDTH)])
    rhs = np.array([(factor * offset**k).deriv(derivative)(x0) for k in range(STENCIL_WIDTH)])
    return np.linalg.solve(matrix, rhs), float(np.linalg.cond(matrix))


def _field_matrix(nodes: np.ndarray, derivative: int, p_north: int, p_south: int) -> tuple[sparse.csr_matrix, float]:
    matrix = _ordinary_matrix(nodes, derivative).tolil()
    conditions = []
    for row in range(3):
        columns = np.arange(STENCIL_WIDTH)
        weights, condition = _factor_weights(nodes[columns], nodes[row], derivative, p_north, True)
        matrix.rows[row] = columns.tolist(); matrix.data[row] = weights.tolist(); conditions.append(condition)
    for row in range(len(nodes) - 3, len(nodes)):
        columns = np.arange(len(nodes) - STENCIL_WIDTH, len(nodes))
        weights, condition = _factor_weights(nodes[columns], nodes[row], derivative, p_south, False)
        matrix.rows[row] = columns.tolist(); matrix.data[row] = weights.tolist(); conditions.append(condition)
    return matrix.tocsr(), max(conditions)


@dataclass(frozen=True)
class SpatialGrid:
    R: np.ndarray
    y: np.ndarray
    RR: np.ndarray
    yy: np.ndarray
    D1R: sparse.csr_matrix
    D2R: sparse.csr_matrix
    D1y_field: sparse.csr_matrix
    D2y_field: sparse.csr_matrix
    D1y_source: sparse.csr_matrix
    D2y_source: sparse.csr_matrix
    dR: float
    dy: float
    p_north: int
    p_south: int
    version: str
    max_factor_condition: float

    def dr(self, values: np.ndarray) -> np.ndarray:
        return self.D1R @ values

    def drr(self, values: np.ndarray) -> np.ndarray:
        return self.D2R @ values

    def dy_field(self, values: np.ndarray) -> np.ndarray:
        return (self.D1y_field @ values.T).T

    def dyy_field(self, values: np.ndarray) -> np.ndarray:
        return (self.D2y_field @ values.T).T

    def dy_source(self, values: np.ndarray) -> np.ndarray:
        return (self.D1y_source @ values.T).T


def build_spatial_grid(chi: float, m: int, n_r: int, n_y: int) -> SpatialGrid:
    if n_r < 7 or n_y < 7 or n_y % 2 == 0:
        raise ContractError("seven-point grid requires n_r,n_y >= 7 and odd n_y for augmented Simpson projection")
    R_horizon = 1 / horizons(1.0, chi)[0]
    R = np.linspace(R_horizon, 0.0, n_r + 1)
    y = np.linspace(-1.0, 1.0, n_y + 2)[1:-1]
    p_north, p_south = pole_powers(m)
    metadata = json.dumps({"stencil": 7, "n_r": n_r, "n_y": n_y, "m": m, "chi": chi}, sort_keys=True)
    RR, yy = np.meshgrid(R, y, indexing="ij")
    D1yf, condition1 = _field_matrix(y, 1, p_north, p_south)
    D2yf, condition2 = _field_matrix(y, 2, p_north, p_south)
    return SpatialGrid(
        R, y, RR, yy, _ordinary_matrix(R, 1), _ordinary_matrix(R, 2),
        D1yf, D2yf,
        _ordinary_matrix(y, 1), _ordinary_matrix(y, 2),
        abs(R[1] - R[0]), y[1] - y[0], p_north, p_south,
        hashlib.sha256(metadata.encode("ascii")).hexdigest(), max(condition1, condition2),
    )


def reconstruct_axis_endpoints(values: np.ndarray, grid: SpatialGrid) -> np.ndarray:
    values = np.asarray(values)
    if values.shape[-1] != grid.y.size:
        raise ContractError("field angular axis does not match grid")
    result = np.empty(values.shape[:-1] + (grid.y.size + 2,), dtype=values.dtype)
    result[..., 1:-1] = values
    for target_index, power, columns in ((0, grid.p_north, np.arange(7)), (-1, grid.p_south, np.arange(grid.y.size - 7, grid.y.size))):
        if power > 0:
            result[..., target_index] = 0
        else:
            endpoint = -1.0 if target_index == 0 else 1.0
            weights = finite_difference_weights(endpoint, grid.y[columns], 0)
            result[..., target_index] = np.tensordot(values[..., columns], weights, axes=(-1, 0))
    return result


def simpson_weights_augmented(grid: SpatialGrid) -> tuple[np.ndarray, np.ndarray]:
    nodes = np.concatenate(([-1.0], grid.y, [1.0]))
    intervals = len(nodes) - 1
    if intervals % 2:
        raise ContractError("composite Simpson requires an even number of intervals")
    weights = np.ones_like(nodes)
    weights[1:-1:2] = 4
    weights[2:-1:2] = 2
    weights *= grid.dy / 3
    return nodes, weights
