import sys
from pathlib import Path
import numpy as np
from scipy.special import pro_cv

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
from environment_cloud import angular_eigenvalue, cloud_211, horizon, mode_flux


def test_angular_spectrum_against_independent_scipy():
    for m in (0, 1, 2):
        for ell in range(m, m+4):
            for c in (0., .2, 1.):
                np.testing.assert_allclose(angular_eigenvalue(ell, m, -c*c),
                                           pro_cv(m, ell, c), atol=1e-12)


def test_cloud_spectrum_and_published_threshold():
    base = cloud_211()
    fine = cloud_211(outer_efolds=45., horizon_offset=1e-6, rtol=2e-11)
    assert abs(base['m2_threshold_r_over_M']-41.66) < .005
    assert abs(base['m2_threshold_r_over_M']-fine['m2_threshold_r_over_M']) < 1e-6
    assert abs(fine['matching_residual']) < 1e-9


def test_charge_bookkeeping_and_horizon_area_law():
    a, mu, omega_p = .87715302759493, .3, .01
    _, omega_c = horizon(a)
    for m in range(-40, 5):  # includes negative-frequency propagating channels
        omega = omega_c+(m-1)*omega_p
        flux = mode_flux(omega, m, omega_c, 1, mu, a, 1+2j, .3-.4j)
        for f in flux.values():
            np.testing.assert_allclose(f['orbital_energy'],
                                       omega_p*f['orbital_angular_momentum'], atol=1e-14)
            np.testing.assert_allclose(f['wave_energy']-omega_c*f['charge'],
                                       f['orbital_energy'], atol=1e-14)
        assert flux['infinity']['wave_energy'] >= 0
        if abs(omega) < mu:
            assert flux['infinity']['charge'] == 0
        f = flux['horizon']
        assert f['wave_energy']-omega_c*m*f['charge'] >= -1e-14
