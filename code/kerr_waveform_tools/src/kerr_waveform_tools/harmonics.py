from __future__ import annotations

import math

import numpy as np

from .contracts import ContractError


def _wigner_d(ell: int, mp: int, m: int, theta) -> np.ndarray:
    theta = np.asarray(theta, dtype=float)
    prefactor = math.sqrt(
        math.factorial(ell + m)
        * math.factorial(ell - m)
        * math.factorial(ell + mp)
        * math.factorial(ell - mp)
    )
    total = np.zeros_like(theta)
    kmin = max(0, m - mp)
    kmax = min(ell + m, ell - mp)
    for k in range(kmin, kmax + 1):
        denominator = (
            math.factorial(ell + m - k)
            * math.factorial(k)
            * math.factorial(mp - m + k)
            * math.factorial(ell - mp - k)
        )
        total += (
            (-1) ** (k - m + mp)
            * prefactor
            / denominator
            * np.cos(theta / 2) ** (2 * ell + m - mp - 2 * k)
            * np.sin(theta / 2) ** (mp - m + 2 * k)
        )
    return total


def spin_weighted_spherical_harmonic(
    ell: int,
    m: int,
    y,
    spin: int = -2,
) -> np.ndarray:
    if ell < abs(spin) or abs(m) > ell:
        raise ContractError("invalid spin-weighted harmonic indices")
    y = np.asarray(y, dtype=float)
    if np.any(np.abs(y) > 1):
        raise ContractError("angular coordinate y must lie in [-1, 1]")
    theta = np.arccos(-y)
    return (
        (-1) ** spin
        * np.sqrt((2 * ell + 1) / (4 * np.pi))
        * _wigner_d(ell, m, -spin, theta)
    )


def projection_row(ell: int, m: int, y, weights) -> np.ndarray:
    y, weights = np.asarray(y, dtype=float), np.asarray(weights, dtype=float)
    if y.shape != weights.shape:
        raise ContractError("quadrature nodes and weights must have equal shape")
    harmonic = spin_weighted_spherical_harmonic(ell, m, y)
    return 2 * np.pi * weights * np.conjugate(harmonic)


def project_mode(scri_field, row) -> complex:
    scri_field, row = np.asarray(scri_field), np.asarray(row)
    if scri_field.shape[-1] != row.shape[0]:
        raise ContractError("last field axis must match projection row")
    return np.sum(scri_field * row, axis=-1)


__all__ = [
    "project_mode",
    "projection_row",
    "spin_weighted_spherical_harmonic",
]
