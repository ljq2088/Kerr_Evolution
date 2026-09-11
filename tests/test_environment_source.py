import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
from environment_source import ThresholdCloud, angular_mode, kerr_metric, connection, project_source


@pytest.fixture(scope='module')
def cloud():
    return ThresholdCloud()


def test_independent_mass_charge_identity_and_quadrature(cloud):
    e, q = cloud.integrals(4801)
    assert abs(e-1) < 1e-7
    assert abs(e-cloud.omega*q) < 1e-7


def test_angular_norm_and_derivative(cloud):
    x, w = np.polynomial.legendre.leggauss(80)
    theta = np.arccos(x)
    s, ds, _ = angular_mode(theta, 1, 1, cloud.c2)
    np.testing.assert_allclose(2*np.pi*np.dot(w, s*s), 1, atol=1e-13)
    sp = angular_mode(theta+1e-5, 1, 1, cloud.c2)[0]
    sm = angular_mode(theta-1e-5, 1, 1, cloud.c2)[0]
    np.testing.assert_allclose((sp-sm)/2e-5, ds, atol=1e-9)


def test_covariant_hessian_trace_and_metric_compatibility(cloud):
    for r in (cloud.rp+.01, 3.5, 20., 41.6, 100.):
        for theta in (.3, 1.2, 2.6):
            phi, hess, inv = cloud.hessian(r, theta)
            np.testing.assert_allclose(np.einsum('ij,ij', inv, hess),
                                       cloud.mu**2*phi, rtol=2e-9, atol=1e-13)
            # A constant conformal test metric is Lorenz; its source = mu² phi.
            np.testing.assert_allclose(cloud.lorenz_source(r, theta, kerr_metric(r,theta,cloud.a)),
                                       cloud.mu**2*phi, rtol=2e-9, atol=1e-13)
            _, gamma = connection(r, theta, cloud.a)
            g = kerr_metric(r, theta, cloud.a)
            dg = (kerr_metric(r+1e-5, theta, cloud.a)-
                  kerr_metric(r-1e-5, theta, cloud.a))/2e-5
            expected = gamma[:, 1, :].T@g+g@gamma[:, 1, :]
            np.testing.assert_allclose(dg, expected, rtol=2e-6, atol=2e-8)


def test_hessian_against_independent_finite_difference(cloud):
    point, step = np.array([.0, 20., 1.2, .0]), 3e-4
    def field(p):
        return (cloud.radial(p[1])[0,0]*angular_mode(p[2],1,1,cloud.c2)[0]
                *np.exp(-1j*cloud.omega*p[0]+1j*p[3]))
    shifts = np.eye(4)*step
    f0 = field(point)
    grad = np.array([(field(point+s)-field(point-s))/(2*step) for s in shifts])
    partial = np.empty((4,4), complex)
    for i in range(4):
        for j in range(4):
            if i == j:
                partial[i,j] = (field(point+shifts[i])-2*f0+field(point-shifts[i]))/step**2
            else:
                u,v = shifts[i], shifts[j]
                partial[i,j] = (field(point+u+v)-field(point+u-v)-field(point-u+v)
                                +field(point-u-v))/(4*step**2)
    _, gamma = connection(20., 1.2, cloud.a)
    numerical = partial-np.einsum('kij,k->ij', gamma, grad)
    np.testing.assert_allclose(cloud.hessian(20.,1.2)[1], numerical, rtol=2e-5, atol=1e-10)


def test_source_projection_includes_sigma(cloud):
    # h=g is a test tensor, not a particle. Source=mu² phi and Sigma couples l=1,3.
    r = np.array([5.,20.])
    callback = lambda r,t: kerr_metric(r,t,cloud.a)
    w, j1 = project_source(cloud,r,20.,1,1,callback)
    _, j3 = project_source(cloud,r,20.,3,1,callback)
    assert w == cloud.omega
    x, weights = np.polynomial.legendre.leggauss(120)
    s1 = angular_mode(np.arccos(x),1,1,cloud.c2)[0]
    s3 = angular_mode(np.arccos(x),3,1,cloud.c2)[0]
    R = cloud.radial(r)[0]
    expected = cloud.mu**2*R*2*np.pi*(r*r*np.dot(weights,s1*s1)
                                    +cloud.a**2*np.dot(weights,x*x*s1*s1))
    np.testing.assert_allclose(j1, expected, rtol=1e-11)
    np.testing.assert_allclose(j3, cloud.mu**2*R*2*np.pi*cloud.a**2*np.dot(weights,x*x*s1*s3),rtol=1e-9)


def test_high_angular_eigenfunctions_share_radial_spectrum_and_normalization():
    from environment_cloud import angular_eigenvalue
    from environment_source import angular_mode
    x, weights = np.polynomial.legendre.leggauss(120)
    theta = np.arccos(x)
    for ell in (21, 24, 32):
        for m in (-1, 1, 5):
            for c2 in (-1., .01, 1.):
                value, derivative, lam = angular_mode(theta, ell, m, c2)
                fine, dfine, _ = angular_mode(theta, ell, m, c2, size=80)
                np.testing.assert_allclose(lam, angular_eigenvalue(ell, m, c2), atol=2e-11, rtol=2e-13)
                np.testing.assert_allclose(2*np.pi*np.dot(weights, value**2), 1., atol=2e-12)
                np.testing.assert_allclose(value, fine, atol=2e-11, rtol=2e-11)
                np.testing.assert_allclose(derivative, dfine, atol=2e-10, rtol=2e-10)
