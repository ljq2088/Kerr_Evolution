"""Massive |211> threshold cloud: independent radial shooting, M=1.

Dyson et al. arXiv:2501.09806v1. This is a background spectral calculation,
not a Lorenz metric reconstruction or a forced environmental waveform.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


def angular_eigenvalue(ell, m, c2, size=20):
    """Eigenvalue of -angular Laplacian - c² cos²(theta), unit sphere norm."""
    ls = np.arange(abs(m), abs(m) + size + 1)
    cosine = np.zeros((size + 1, size + 1))
    for j, l in enumerate(ls[:-1]):
        cosine[j, j+1] = cosine[j+1, j] = np.sqrt(
            ((l+1)**2-m*m)/((2*l+1)*(2*l+3)))
    matrix = np.diag(ls[:size]*(ls[:size]+1.)) - c2*(cosine@cosine)[:size, :size]
    return float(np.linalg.eigvalsh(matrix)[ell-abs(m)])


def horizon(a):
    if not 0 < a < 1:
        raise ValueError('Require 0 < a/M < 1')
    rp = 1 + np.sqrt(1-a*a)
    return rp, a/(2*rp)


def radial_coefficients(r, a, mu, omega, m, lam):
    """Delta, Delta', V for (Delta R')' + V R = J, M=1."""
    delta = r*r-2*r+a*a
    K = (r*r+a*a)*omega-a*m
    return delta, 2*(r-1), K*K/delta-mu*mu*r*r-a*a*omega*omega+2*a*m*omega-lam


def matching_residual(a, mu=.3, outer_efolds=30., horizon_offset=1e-5,
                      rtol=2e-10):
    """Match horizon-regular and infinity-decaying solutions at threshold.

    At omega=m Omega_H, K(r+) vanishes; the regular Frobenius branch has
    R'(r+)=-V(r+)R(r+)/Delta'(r+). Infinity uses the first Coulomb power.
    Finite boundary errors must be checked by changing both cutoffs.
    """
    rp, omega = horizon(a)  # m_b=1
    if omega >= mu:
        raise ValueError('Cloud must decay: omega < mu')
    lam = angular_eigenvalue(1, 1, a*a*(omega*omega-mu*mu))
    k = np.sqrt(mu*mu-omega*omega)
    def rhs(r, state):
        delta, dp, potential = radial_coefficients(r, a, mu, omega, 1, lam)
        return [state[1], -(dp*state[1]+potential*state[0])/delta]
    vh = -mu*mu*rp*rp-a*a*omega*omega+2*a*omega-lam
    slope = -vh/(2*(rp-1))
    match = 2/mu**2
    rmax = outer_efolds/k
    beta = (2*omega*omega-mu*mu)/k-1
    left = solve_ivp(rhs, (rp+horizon_offset, match),
                     [1+slope*horizon_offset, slope], rtol=rtol, atol=rtol*1e-2)
    right = solve_ivp(rhs, (rmax, match), [1., -k+beta/rmax],
                      rtol=rtol, atol=rtol*1e-2)
    if not left.success or not right.success:
        raise RuntimeError('Radial integration failed')
    u, v = left.y[:, -1], right.y[:, -1]
    return float((u[0]*v[1]-v[0]*u[1])/np.linalg.norm(u)/np.linalg.norm(v))


def cloud_211(outer_efolds=30., horizon_offset=1e-5, rtol=2e-10):
    """Find the alpha=.3 |211> threshold branch; bracket is deliberately fixed."""
    kwargs = dict(outer_efolds=outer_efolds, horizon_offset=horizon_offset, rtol=rtol)
    a = brentq(lambda x: matching_residual(x, **kwargs), .87, .88, xtol=1e-12)
    rp, omega = horizon(a)
    mu = .3
    threshold = (1/(mu-omega)-a)**(2/3)
    return dict(alpha=mu, a_over_M=a, M_omega_c=omega,
                r_plus_over_M=rp, m2_threshold_r_over_M=threshold,
                matching_residual=matching_residual(a, **kwargs))


def mode_flux(omega, m, omega_c, mb, mu, a, z_inf, z_h):
    """C=1, angular norm=1; omit common epsilon² q². Inward horizon positive.

    Return ordinary wave and effective orbital fluxes separately. Massive
    infinity flux is a stationary large-radius limit, not massless scri data.
    """
    rp, oh = horizon(a)
    ni = 2*np.sign(omega)*np.sqrt(max(0., omega*omega-mu*mu))*abs(z_inf)**2
    nh = 2*(rp*rp+a*a)*(omega-m*oh)*abs(z_h)**2
    return {name: dict(charge=float(n), wave_energy=float(omega*n),
                       orbital_energy=float((omega-omega_c)*n),
                       orbital_angular_momentum=float((m-mb)*n))
            for name, n in [('infinity', ni), ('horizon', nh)]}


if __name__ == '__main__':
    baseline = cloud_211()
    refined = cloud_211(outer_efolds=45., horizon_offset=1e-6, rtol=2e-11)
    report = dict(method='Independent two-sided shooting; not Leaver',
                  baseline=baseline, refined=refined,
                  paper_v1_m2_threshold=41.66,
                  threshold_cutoff_change=abs(baseline['m2_threshold_r_over_M']-
                                              refined['m2_threshold_r_over_M']))
    destination = Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    destination.mkdir(exist_ok=True)
    (destination/'cloud_211.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))
