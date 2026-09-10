import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
from environment_radial import RadialGreen
from environment_cloud import radial_coefficients


@pytest.mark.parametrize('omega', [.4, -.4, .2])
def test_green_manufactured_compact_solution(omega):
    # Nonzero exact interior field, zero at both physical asymptotes.
    solver = RadialGreen(.7, .3, omega, 2, 2, rmax=120.)
    errors = []
    for n in (401, 801):
        # Place the piecewise-smooth source boundaries on mesh points.
        r = np.r_[np.linspace(solver.rmin, 10., 201)[:-1],
                  np.linspace(10., 30., n), np.linspace(30., solver.rmax, 201)[1:]]
        x = (r-10)/20
        inside = (x > 0) & (x < 1)
        p = np.polynomial.Polynomial([0, 0, 0, 0, 1, -4, 6, -4, 1])*256
        f = np.where(inside, p(x), 0.)
        df = np.where(inside, p.deriv()(x)/20, 0.)
        ddf = np.where(inside, p.deriv(2)(x)/400, 0.)
        d, dp, v = radial_coefficients(r, .7, .3, omega, 2, solver.lam)
        result = solver.solve(r, d*ddf+dp*df+v*f)
        errors.append(np.max(np.abs(result['field']-f)))
        assert result['wronskian_relative_spread'] < 1e-6
    assert errors[-1] < 3e-6
    assert errors[-1] < errors[0]/5


def test_threshold_rejected():
    with pytest.raises(ValueError, match='threshold'):
        RadialGreen(.7, .3, .3, 2, 2)
