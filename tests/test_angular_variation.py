import numpy as np
import pytest
from environment_source import angular_mode
from environment_angular_variation import angular_mode_mass_derivative


@pytest.mark.parametrize('ell,m,a,omega',[(1,1,.877,.03),(18,1,.877,.03),(24,2,.7,.2),(18,-1,.877,.03)])
def test_mass_derivative_matches_separate_eigensolves_and_norm(ell,m,a,omega):
    x,w=np.polynomial.legendre.leggauss(90);theta=np.arccos(x)
    s,ds,lam,dm,ddm,dlam=angular_mode_mass_derivative(theta,ell,m,a,omega)
    base=angular_mode(theta,ell,m,a*a*omega*omega)
    np.testing.assert_allclose(s,base[0],rtol=0,atol=2e-13)
    np.testing.assert_allclose(ds,base[1],rtol=0,atol=2e-12)
    assert abs(lam-base[2])<1e-12
    h=.01
    def shifted(t):return angular_mode(theta,ell,m,a*a*(omega*omega-t))
    plus,minus,halfplus,halfminus=[shifted(t) for t in (h,-h,h/2,-h/2)]
    for index,actual in enumerate((dm,ddm,dlam)):
        coarse=(plus[index]-minus[index])/(2*h)
        fine=(halfplus[index]-halfminus[index])/h
        np.testing.assert_allclose(actual,(4*fine-coarse)/3,rtol=2e-6,atol=2e-9)
    assert abs(2*np.pi*np.dot(w,s*s)-1)<2e-12
    assert abs(2*np.pi*np.dot(w,s*dm))<2e-13
    assert abs(dlam-a*a*2*np.pi*np.dot(w,x*x*s*s))<2e-12


def test_spherical_limit_and_axis_guard():
    result=angular_mode_mass_derivative(np.array([.4,1.2,2.1]),3,-1,0.,.2)
    for value in result[3:]:assert np.all(value==0)
    with pytest.raises(ValueError,match='Polar derivative'):
        angular_mode_mass_derivative(0.,1,1,.8,.1)
