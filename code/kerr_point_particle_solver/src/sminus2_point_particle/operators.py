from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np

from .errors import ContractError, OpenBoundaryError


def uniform_nodes(start: float, stop: float, count: int) -> np.ndarray:
    if count < 3 or not start < stop:
        raise ContractError("uniform grid requires count >= 3 and start < stop")
    return np.linspace(start, stop, count)


def open_uniform_nodes(start: float, stop: float, count: int) -> np.ndarray:
    if count < 3 or not start < stop:
        raise ContractError("open grid requires count >= 3 and start < stop")
    return np.linspace(start, stop, count + 2)[1:-1]


def finite_difference_weights(x0: float, nodes: np.ndarray, derivative: int) -> np.ndarray:
    """Fornberg weights formed on a locally scaled stencil."""
    nodes = np.asarray(nodes, dtype=float)
    n = len(nodes)
    if derivative < 0 or derivative >= n or len(np.unique(nodes)) != n:
        raise ContractError("invalid derivative order or repeated stencil nodes")
    scale = np.max(np.abs(nodes - x0))
    if scale == 0:
        raise ContractError("finite-difference stencil has zero scale")
    x = (nodes - x0) / scale
    weights = np.zeros((n, derivative + 1))
    weights[0, 0] = 1.0
    c1 = 1.0
    c4 = x[0]
    for i in range(1, n):
        mn = min(i, derivative)
        c2 = 1.0
        c5 = c4
        c4 = x[i]
        for j in range(i):
            c3 = x[i] - x[j]
            c2 *= c3
            if j == i - 1:
                for k in range(mn, 0, -1):
                    weights[i, k] = c1 * (k * weights[i - 1, k - 1] - c5 * weights[i - 1, k]) / c2
                weights[i, 0] = -c1 * c5 * weights[i - 1, 0] / c2
            for k in range(mn, 0, -1):
                weights[j, k] = (c4 * weights[j, k] - k * weights[j, k - 1]) / c3
            weights[j, 0] = c4 * weights[j, 0] / c3
        c1 = c2
    return weights[:, derivative] / scale**derivative


@dataclass(frozen=True)
class InteriorDerivative:
    nodes: np.ndarray
    matrix: np.ndarray
    valid_rows: np.ndarray
    derivative: int

    def apply(self, values: np.ndarray, axis: int = 0) -> np.ndarray:
        values = np.asarray(values)
        if values.shape[axis] != len(self.nodes):
            raise ContractError("derivative axis length does not match nodes")
        moved = np.moveaxis(values, axis, 0)
        out = np.full(moved.shape, np.nan, dtype=np.result_type(values, float))
        for row in self.valid_rows:
            columns = np.flatnonzero(self.matrix[row])
            out[row] = np.tensordot(self.matrix[row, columns], moved[columns], axes=(0, 0))
        return np.moveaxis(out, 0, axis)


def interior_derivative(nodes: np.ndarray, derivative: int, stencil_size: int = 7) -> InteriorDerivative:
    nodes = np.asarray(nodes, dtype=float)
    if stencil_size % 2 == 0 or stencil_size <= derivative or stencil_size > len(nodes):
        raise ContractError("use an odd stencil larger than derivative order")
    half = stencil_size // 2
    matrix = np.full((len(nodes), len(nodes)), np.nan)
    rows = np.arange(half, len(nodes) - half)
    for i in rows:
        sl = slice(i - half, i + half + 1)
        matrix[i] = 0.0
        matrix[i, sl] = finite_difference_weights(nodes[i], nodes[sl], derivative)
    return InteriorDerivative(nodes, matrix, rows, derivative)


@dataclass(frozen=True)
class BoundaryRowExtension:
    def require_production(self) -> None:
        raise OpenBoundaryError("C0 contains no production boundary rows; V1.2 must supply and review actual closures")
