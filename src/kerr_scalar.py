"""Massless scalar on Kerr: horizon-to-scri hyperboloidal evolution.

See docs/kerr_hyperboloidal_zh.md for coordinates and all coefficients.
"""
from dataclasses import dataclass, asdict
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.integrate import solve_ivp


def chebyshev(n):
    """Ascending Lobatto nodes on [0,1] and their differentiation matrix."""
    x = -np.cos(np.pi * np.arange(n + 1) / n)
    w = (-1.) ** np.arange(n + 1)
    w[[0, -1]] *= .5
    dx = x[:, None] - x[None, :]
    np.fill_diagonal(dx, 1.)
    d = w[None, :] / w[:, None] / dx
    np.fill_diagonal(d, 0.)
    np.fill_diagonal(d, -d.sum(axis=1))
    return (x + 1) / 2, 2 * d


def sin_squared_matrix(m, lmax):
    """Exact Galerkin projection, including the lmax+1 intermediate mode."""
    ls = np.arange(abs(m), lmax + 2)
    x = np.zeros((len(ls), len(ls)))
    for j, l in enumerate(ls[:-1]):
        x[j, j + 1] = x[j + 1, j] = np.sqrt(((l + 1)**2 - m*m) / (4*(l + 1)**2 - 1))
    return (np.eye(len(ls)) - x @ x)[:-1, :-1]


@dataclass
class Config:
    mass: float = 1.
    spin: float = .7  # a/M
    m: int = 0
    lmax: int = 4
    initial_l: int = 2
    n: int = 64
    tmax: float = 80.
    samples: int = 401
    center: float = .45
    width: float = .12
    rtol: float = 1.e-9
    atol: float = 1.e-11

    def validate(self):
        vals = (self.mass, self.spin, self.tmax, self.center, self.width, self.rtol, self.atol)
        if not np.all(np.isfinite(vals)):
            raise ValueError('parameters must be finite')
        if self.mass <= 0 or abs(self.spin) >= 1:
            raise ValueError('require M > 0 and |a/M| < 1')
        if self.n < 8 or self.samples < 2 or self.tmax <= 0:
            raise ValueError('require n >= 8, samples >= 2, tmax > 0')
        if not abs(self.m) <= self.initial_l <= self.lmax:
            raise ValueError('require |m| <= initial_l <= lmax')
        if not 0 < self.center < 1 or self.width <= 0 or min(self.rtol, self.atol) <= 0:
            raise ValueError('invalid initial profile or tolerance')


class KerrScalar:
    def __init__(self, config):
        config.validate()
        self.config = config
        M, a = config.mass, config.mass * config.spin
        L = M + np.sqrt(M*M - a*a)
        self.horizon = L
        self.s, self.d = chebyshev(config.n)
        self.d2 = self.d @ self.d
        self.ls = np.arange(abs(config.m), config.lmax + 1)
        self.shape = (config.n + 1, len(self.ls))
        s = self.s
        H = 2 + 4*M*s/L
        ar = H * (-8*M*M + 4*M*a*a*s/L)
        A = ar[:, None, None] * np.eye(len(self.ls)) + a*a*sin_squared_matrix(config.m, config.lmax)
        if np.max(np.linalg.eigvalsh(A)) >= 0:
            raise ValueError('time slices must be spacelike')
        self.ainv = np.linalg.inv(A)
        self.C = (2*L + (2*a*a-16*M*M)*s*s/L + 8*M*a*a*s**3/L**2)[:, None]
        self.D = (s*s*(1-2*M*s/L+a*a*s*s/L**2))[:, None]
        self.E = (((2*a*a-16*M*M)*s/L+12*M*a*a*s*s/L**2) + 2j*a*config.m*(1-H))[:, None]
        self.F = (2*s-6*M*s*s/L+4*a*a*s**3/L**2 - 2j*a*config.m*s*s/L)[:, None]
        self.V = (-2*M*s/L+2*a*a*s*s/L**2-2j*a*config.m*s/L)[:, None] - self.ls*(self.ls+1)

    def rhs(self, time, state):
        u, p = state.reshape((2,) + self.shape)
        rest = self.C*(self.d @ p) + self.D*(self.d2 @ u) + self.E*p + self.F*(self.d @ u) + self.V*u
        acceleration = -np.einsum('ijk,ik->ij', self.ainv, rest)
        return np.stack((p, acceleration)).ravel()

    def initial_state(self):
        state = np.zeros((2,) + self.shape, dtype=complex)
        s, c = self.s, self.config
        # Analytic profile with exact zero value and first derivative at endpoints.
        profile = (s*(1-s)/(c.center*(1-c.center)))**2 * np.exp(-((s-c.center)/c.width)**2)
        state[0, :, c.initial_l-abs(c.m)] = profile
        return state.ravel()

    def evolve(self):
        c = self.config
        times = np.linspace(0, c.tmax, c.samples)
        sol = solve_ivp(self.rhs, (0, c.tmax), self.initial_state(), method='DOP853',
                        t_eval=times, rtol=c.rtol, atol=c.atol, max_step=.5*c.mass)
        if not sol.success or not np.all(np.isfinite(sol.y)):
            raise RuntimeError(sol.message)
        fields = sol.y.T.reshape((-1, 2) + self.shape)
        return times, fields, sol.nfev


def save_run(config, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    model = KerrScalar(config)
    times, fields, nfev = model.evolve()
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    u = fields[:, 0]
    np.savez_compressed(out/'evolution.npz', tau=times, sigma=model.s, ell=model.ls,
                        u=u, p=fields[:, 1], config=json.dumps(asdict(config)))
    meta = dict(config=asdict(config), horizon=model.horizon, nfev=nfev,
                max_abs_u=float(np.max(np.abs(u))), final_max_abs_u=float(np.max(np.abs(u[-1]))))
    (out/'run.json').write_text(json.dumps(meta, indent=2)+'\n')
    columns = [times]
    names = ['tau']
    for endpoint, index in [('scri', 0), ('horizon', -1)]:
        for j, l in enumerate(model.ls):
            columns.extend([u[:, index, j].real, u[:, index, j].imag])
            names.extend([f'{endpoint}_l{l}_real', f'{endpoint}_l{l}_imag'])
    np.savetxt(out/'waveforms.csv', np.column_stack(columns), delimiter=',', header=','.join(names), comments='')
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for ax, index, title in zip(axes, [0, -1], ['Future null infinity', 'Future event horizon']):
        for j, l in enumerate(model.ls):
            if np.max(np.abs(u[:, index, j])) > 1e-13:
                ax.plot(times/config.mass, u[:, index, j].real, label=f'l={l}, Re(u)')
                if config.m:
                    ax.plot(times/config.mass, u[:, index, j].imag, '--', label=f'l={l}, Im(u)')
        ax.set(xlabel='tau / M', ylabel='u = r Phi (mode coefficient)', title=title)
        ax.legend()
        ax.grid(alpha=.25)
    fig.suptitle(f'Kerr scalar: a/M={config.spin}, m={config.m}')
    fig.savefig(out/'waveforms.png', dpi=180)
    plt.close(fig)
    print(json.dumps(meta, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name, field in Config.__dataclass_fields__.items():
        parser.add_argument('--'+name.replace('_', '-'), type=type(field.default), default=field.default)
    parser.add_argument('--output', default='outputs/kerr_scalar')
    args = vars(parser.parse_args())
    output = args.pop('output')
    save_run(Config(**args), output)
