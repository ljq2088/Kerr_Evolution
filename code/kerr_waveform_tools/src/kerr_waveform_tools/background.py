from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .contracts import ContractError


def delta(r, M, a):
    return np.asarray(r) ** 2 - 2 * M * np.asarray(r) + a**2


def sigma(r, y, a):
    return np.asarray(r) ** 2 + a**2 * np.asarray(y) ** 2


def horizons(M: float, a: float) -> tuple[float, float]:
    if M <= 0 or abs(a) >= M:
        raise ContractError("C0 uses nonextremal Kerr: M > 0 and |a| < M")
    root = np.sqrt(M * M - a * a)
    return M + root, M - root


def compactified_radius(r, L):
    r = np.asarray(r)
    if np.any(r <= 0) or L <= 0:
        raise ContractError("r and L must be positive")
    return L**2 / r


def boyer_lindquist_radius(R, L):
    R = np.asarray(R)
    if np.any(R <= 0) or L <= 0:
        raise ContractError("R and L must be positive; R=0 is handled only by SCRI interfaces")
    return L**2 / R


def compact_delta(R, M, a, L):
    R = np.asarray(R)
    return L**4 - 2 * M * L**2 * R + a**2 * R**2


def compact_sigma(R, y, a, L):
    return L**4 + a**2 * np.asarray(R) ** 2 * np.asarray(y) ** 2


def isco_quantities(chi: float, M: float = 1.0) -> tuple[float, float, float, float]:
    if not 0 <= chi < 1 or M <= 0:
        raise ContractError("baseline ISCO formula requires M > 0 and 0 <= chi < 1")
    z1 = 1 + (1 - chi**2) ** (1 / 3) * ((1 + chi) ** (1 / 3) + (1 - chi) ** (1 / 3))
    z2 = np.sqrt(3 * chi**2 + z1**2)
    rt = 3 + z2 - np.sqrt((3 - z1) * (3 + z1 + 2 * z2))
    v = rt ** -0.5
    den = np.sqrt(1 - 3 * v**2 + 2 * chi * v**3)
    energy = (1 - 2 * v**2 + chi * v**3) / den
    lz = M * (1 - 2 * chi * v**3 + chi**2 * v**4) / (v * den)
    omega = 1 / (M * (rt**1.5 + chi))
    return M * rt, float(energy), float(lz), float(omega)


@dataclass(frozen=True)
class RadialState:
    r: float
    potential: float
    potential_sqrt: float
    K: float
    u_t_bl: float
    u_r_bl: float
    u_phi_bl: float
    u_T: float
    u_R: float
    u_Phi: float
    u_T_r: float
    u_T_rr: float
    u_R_r: float
    u_R_rr: float
    u_Phi_r: float
    u_Phi_rr: float
    K_r: float
    K_rr: float


def radial_state(r: float, M: float, a: float, L: float, energy: float, angular_momentum: float, *, allow_horizon_interior: bool = False) -> RadialState:
    if L <= 0:
        raise ContractError("compactification scale L must be positive")
    r_plus = horizons(M, a)[0]
    if r <= r_plus and not allow_horizon_interior:
        raise ContractError("C0 radial jets are restricted to r > r_+")
    B = angular_momentum - a * energy
    P = energy * (r * r + a * a) - a * angular_momentum
    C = r * r + B * B
    de = float(delta(r, M, a))
    rad = P * P - de * C
    if rad <= 0:
        raise ContractError("A5 radial-jet branch requires strictly positive radial potential")
    V = np.sqrt(rad)
    H = P + V
    if H == 0:
        raise ContractError("A5 regular K branch requires P + sqrt(R) != 0")
    P1, P2 = 2 * energy * r, 2 * energy
    d1, d2 = 2 * (r - M), 2.0
    C1 = C2 = 2 * r
    C2 = 2.0
    rad1 = 2 * P * P1 - d1 * C - de * C1
    rad2 = 2 * (P1 * P1 + P * P2) - d2 * C - 2 * d1 * C1 - de * C2
    V1 = rad1 / (2 * V)
    V2 = rad2 / (2 * V) - rad1**2 / (4 * V**3)
    H1, H2 = P1 + V1, P2 + V2
    K = C / H
    K1 = C1 / H - C * H1 / H**2
    K2 = C2 / H - 2 * C1 * H1 / H**2 - C * H2 / H**2 + 2 * C * H1**2 / H**3
    eta = 2 + 4 * M / r
    eta1, eta2 = -4 * M / r**2, 8 * M / r**3
    FT = a * B + (r * r + a * a) * K + eta * V
    FT1 = 2 * r * K + (r * r + a * a) * K1 + eta1 * V + eta * V1
    FT2 = 2 * K + 4 * r * K1 + (r * r + a * a) * K2 + eta2 * V + 2 * eta1 * V1 + eta * V2
    FP, FP1, FP2 = B + a * K, a * K1, a * K2
    quotient = lambda F, F1, F2: (F / r**2, F1 / r**2 - 2 * F / r**3, F2 / r**2 - 4 * F1 / r**3 + 6 * F / r**4)
    uT, uT1, uT2 = quotient(FT, FT1, FT2)
    uP, uP1, uP2 = quotient(FP, FP1, FP2)
    uR = L**2 * V / r**4
    uR1 = L**2 * (V1 / r**4 - 4 * V / r**5)
    uR2 = L**2 * (V2 / r**4 - 8 * V1 / r**5 + 20 * V / r**6)
    ur = -V / r**2
    if de > 0:
        ut = (a * B + (r * r + a * a) * P / de) / r**2
        uphi = (B + a * P / de) / r**2
    else:
        ut = uphi = np.nan
    return RadialState(r, rad, V, K, ut, ur, uphi, uT, uR, uP, uT1, uT2, uR1, uR2, uP1, uP2, K1, K2)


def bl_timelike_norm(state: RadialState, M: float, a: float) -> float:
    r = state.r
    de = float(delta(r, M, a))
    if de <= 0:
        raise ContractError("BL timelike norm is not evaluated at or inside the future horizon")
    gtt = -(1 - 2 * M / r)
    gtphi = -2 * M * a / r
    grr = r * r / de
    gphiphi = r * r + a * a + 2 * M * a * a / r
    return gtt * state.u_t_bl**2 + 2 * gtphi * state.u_t_bl * state.u_phi_bl + grr * state.u_r_bl**2 + gphiphi * state.u_phi_bl**2
