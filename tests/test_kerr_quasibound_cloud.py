"""Focused invariants for the bounded complex Kerr cloud implementation."""
import numpy as np
import pytest
from environment_kerr_quasibound_cloud import GeneralKerrQuasiboundCloud,ComplexScalarAngular
from environment_radial import horizon_series

@pytest.fixture(scope="module")
def cloud():
    return GeneralKerrQuasiboundCloud()

def test_complex_frequency_and_KS_mass_balance(cloud):
    assert cloud.omega.imag>0 and cloud.omega==cloud.spectral_omega
    assert abs(cloud.ks_integrals(12001)[0]-1)<1e-8
    balance=cloud.horizon_balance(12001)
    assert balance["horizon_energy"]<0
    assert balance["energy_balance_relative"]<1e-8
    assert balance["charge_balance_relative"]<1e-8

def test_ingoing_BL_phase_matches_independent_horizon_series(cloud):
    expected=np.array(horizon_series(cloud.rmin,cloud.a,cloud.mu,
       cloud.spectral_omega,cloud.m,cloud.lam,order=6))
    actual=cloud.radial(cloud.rmin)[:,0]/cloud.amplitude
    assert np.linalg.norm(actual-expected)/np.linalg.norm(expected)<1e-10

@pytest.mark.parametrize("invalid",[np.nan,np.inf,-np.inf,0.,1e10])
def test_domain_guard_precedes_dense_solution(cloud,invalid,monkeypatch):
    def forbidden(*args):raise AssertionError("Dense output must not be evaluated")
    monkeypatch.setattr(cloud.left,"sol",forbidden)
    monkeypatch.setattr(cloud.right,"sol",forbidden)
    with pytest.raises(ValueError,match="outside solved"):
        cloud.radial(invalid)

def test_complex_angular_uses_bilinear_dual():
    x,w=np.polynomial.legendre.leggauss(64);theta=np.arccos(x)
    one=ComplexScalarAngular(1,1,.1+.2j);three=ComplexScalarAngular(3,1,.1+.2j)
    s1=one.evaluate(theta)[0];s3=three.evaluate(theta)[0]
    assert abs(2*np.pi*np.dot(w,one.dual(theta)*s1)-1)<1e-12
    assert abs(2*np.pi*np.dot(w,one.dual(theta)*s3))<1e-12
    assert abs(2*np.pi*np.dot(w,s1.conjugate()*s3))>1e-3

def test_freezing_preserves_spatial_profile_and_exposes_KG_defect(cloud):
    frozen=GeneralKerrQuasiboundCloud(freeze_growth=True)
    np.testing.assert_allclose(frozen.radial([3,20,100]),cloud.radial([3,20,100]),rtol=0,atol=0)
    assert isinstance(frozen.omega,float) and frozen.spectral_omega.imag>0
    phi,H,inv=frozen.hessian(20.,1.1)
    residual=np.einsum("ij,ij->",inv,H)-frozen.mu**2*phi
    w=frozen.spectral_omega;wr=frozen.omega
    expected=(inv[0,0]*(w*w-wr*wr)+2*inv[0,3]*frozen.m*(wr-w))*phi
    assert abs(residual)/(frozen.mu**2*abs(phi))>1e-9
    assert abs(residual-expected)/(frozen.mu**2*abs(phi))<1e-12

def test_jet_spatial_derivatives_and_KG_equation(cloud):
    g,jet=cloud.jet(20.,1.1,order=4,return_geometry=True)
    R,Rp,Rpp=cloud.radial_state(20.)[:,0]
    S=cloud.angular(1.1)[0]
    np.testing.assert_allclose([jet.value,jet.derivative_value(1,0),jet.derivative_value(2,0)],
      [R*S,Rp*S,Rpp*S],rtol=1e-11,atol=1e-14)
    assert abs(g.scalar_wave(jet).value-cloud.mu**2*jet.value)/(cloud.mu**2*abs(jet.value))<1e-11
