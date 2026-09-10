import sys
from pathlib import Path
import unittest
import numpy as np
from scipy.special import lpmv, gammaln
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from kerr_scalar import Config, KerrScalar, chebyshev, sin_squared_matrix


class KerrTests(unittest.TestCase):
    def test_radial_polynomials(self):
        x, d = chebyshev(24)
        np.testing.assert_allclose(d @ x**5, 5*x**4, atol=2e-12)
        np.testing.assert_allclose(d @ d @ x**5, 20*x**3, atol=2e-9)

    def test_angular_projection_by_quadrature(self):
        x, w = np.polynomial.legendre.leggauss(40)
        for m in (0, 1, 2, 3):
            ls = np.arange(m, 7)
            basis = np.array([lpmv(m, l, x)*np.sqrt((2*l+1)/2*np.exp(gammaln(l-m+1)-gammaln(l+m+1))) for l in ls])
            expected = (basis*(w*(1-x*x))) @ basis.T
            np.testing.assert_allclose(sin_squared_matrix(m, 6), expected, atol=2e-14)

    def test_ripley_equation_11(self):
        # Independent transcription of the published s=0 equation in rho=1/r.
        # Our equation has opposite overall sign and sigma=r_plus*rho.
        for spin in (0., .7, -.8, .99):
            c = Config(mass=1.7, spin=spin, m=2, lmax=5, n=16)
            model = KerrScalar(c)
            M, a, L = c.mass, c.mass*c.spin, model.horizon
            rho = model.s/L
            S = sin_squared_matrix(c.m, c.lmax)
            A = ((16*M*M+8*M*(4*M*M-a*a)*rho-16*a*a*M*M*rho*rho)[:, None, None]*np.eye(len(model.ls)) - a*a*S)
            np.testing.assert_allclose(np.linalg.inv(model.ainv), -A, atol=1e-12)
            np.testing.assert_allclose(model.C[:, 0], 2*L*(1+(a*a-8*M*M)*rho*rho+4*a*a*M*rho**3), atol=1e-13)
            np.testing.assert_allclose(model.D[:, 0], L*L*rho*rho*(1-2*M*rho+a*a*rho*rho), atol=1e-14)
            np.testing.assert_allclose(model.E[:, 0], -2*((8*M*M-a*a)*rho-6*M*a*a*rho*rho+1j*a*c.m*(1+4*M*rho)), atol=1e-13)
            np.testing.assert_allclose(model.F[:, 0], -2*L*rho*(-1+3*M*rho-2*a*a*rho*rho+1j*a*c.m*rho), atol=1e-13)
            expected_v = -2*(M-a*a*rho+1j*a*c.m)*rho
            np.testing.assert_allclose(model.V + model.ls*(model.ls+1), np.broadcast_to(expected_v[:, None], model.shape), atol=1e-14)

    def test_characteristic_endpoints(self):
        for spin in (0., .7, .99, -.99):
            model = KerrScalar(Config(spin=spin))
            A = np.linalg.inv(model.ainv)
            self.assertTrue(np.all(np.linalg.eigvalsh(A) < 0))
            np.testing.assert_allclose(model.D[[0, -1]], 0, atol=1e-15)
            # Principal equation A u_tt + C u_ts + D u_ss = 0:
            # at endpoints one speed is zero, the other is C/A.
            self.assertGreater(model.C[0, 0], 0)  # C/A < 0: leaves sigma >= 0
            self.assertLess(model.C[-1, 0], 0)  # C/A > 0: leaves sigma <= 1

    def test_schwarzschild_mode_decoupling(self):
        model = KerrScalar(Config(spin=0, n=24, tmax=8, samples=4))
        _, fields, _ = model.evolve()
        np.testing.assert_array_equal(fields[:, :, :, [0, 1, 3, 4]], 0)

    def test_m_conjugation_and_kerr_mixing(self):
        plus = KerrScalar(Config(m=1, initial_l=1, n=24, tmax=8, samples=4))
        minus = KerrScalar(Config(m=-1, initial_l=1, n=24, tmax=8, samples=4))
        _, up, _ = plus.evolve()
        _, um, _ = minus.evolve()
        np.testing.assert_allclose(up.conj(), um, atol=1e-12)
        self.assertGreater(np.max(np.abs(up[:, 0, :, 2])), 1e-5)
        np.testing.assert_array_equal(up[:, :, :, [1, 3]], 0)

    def test_mass_rescaling(self):
        one = KerrScalar(Config(mass=1., n=16, tmax=4., samples=5))
        two = KerrScalar(Config(mass=2., n=16, tmax=8., samples=5))
        _, a, _ = one.evolve()
        _, b, _ = two.evolve()
        np.testing.assert_allclose(a[:, 0], b[:, 0], atol=2e-9)
        np.testing.assert_allclose(a[:, 1], 2*b[:, 1], atol=2e-9)

    def test_invalid_inputs(self):
        for kwargs in ({'spin':1.}, {'spin':float('nan')}, {'mass':0.}, {'width':0.}, {'n':4}, {'m':3, 'initial_l':2}):
            with self.assertRaises(ValueError):
                KerrScalar(Config(**kwargs))


if __name__ == '__main__':
    unittest.main()
