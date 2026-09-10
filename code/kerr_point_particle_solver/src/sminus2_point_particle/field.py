from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .errors import ContractError


@dataclass(frozen=True)
class FieldCoefficients:
    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    D: np.ndarray
    E_R: np.ndarray
    F: np.ndarray


def coefficients(R, y, M: float, a: float, L: float, m: int) -> FieldCoefficients:
    R, y = np.broadcast_arrays(np.asarray(R, dtype=float), np.asarray(y, dtype=float))
    d = 1 - y**2
    A = 8 * M * (2 * M - a**2 * R / L**2) * (1 + 2 * M * R / L**2) - a**2 * d
    B = -2 * (L**2 - (8 * M**2 - a**2) * R**2 / L**2 + 4 * a**2 * M * R**3 / L**4)
    C = -(L**2 - 2 * M * R + a**2 * R**2 / L**2) * R**2 / L**2
    D = 2j * a * m * (1 + 4 * M * R / L**2) + 2 * (2 * M * (2 - 3 * a**2 * R**2 / L**4) - a**2 * R / L**2 + 2j * a * y)
    E = 2j * a * m * R**2 / L**2 + 2 * R * (1 + M * R / L**2 - 2 * a**2 * R**2 / L**4)
    F = 2j * a * m * R / L**2 - 2 * M * R / L**2 - 2 * a**2 * R**2 / L**4
    return FieldCoefficients(A, B, C, D, E, F)


def angular_operator(psi, psi_y, psi_yy, y, m: int):
    y = np.asarray(y, dtype=float)
    if np.any(np.abs(y) >= 1):
        raise ContractError("raw A6 angular formula is interior-only; axis closure is open")
    d = 1 - y**2
    return d * psi_yy - 2 * y * psi_y - (m + 2 * y) ** 2 * psi / d - 2 * psi


def p_from_pi(psi, pi, q, coeff: FieldCoefficients):
    return coeff.A * pi + coeff.B * q + coeff.D * psi


def pi_from_p(psi, p, q, coeff: FieldCoefficients):
    return (p - coeff.B * q - coeff.D * psi) / coeff.A

