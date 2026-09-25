import mpmath as mp
import numpy as np
import pytest
from li_leaver_cloud import LiLeaverCloud
from environment_kerr_quasibound_cloud import GeneralKerrQuasiboundCloud

def test_dipole_profile_against_independent_radial_shooting():
    leaver=LiLeaverCloud(terms=300)
    shooting=GeneralKerrQuasiboundCloud()
    assert abs(complex(leaver.omega)-shooting.spectral_omega)<1e-13
    for r in (2.,20.,100.):
        R,Rp=shooting.radial(r)[:,0]
        ratio=complex(leaver.radial_mp(r,1)/leaver.radial_mp(r))
        assert abs(ratio-Rp/R)<2e-9

@pytest.mark.parametrize('ell',[1,2])
def test_finite_series_satisfies_independent_radial_equation(ell):
    c=LiLeaverCloud(ell=ell,m=ell,terms=150)
    assert c.omega.real<c.mu and c.omega.imag>0
    assert abs(c.data['residual'])<mp.mpf('1e-45')
    for r in (3.,20.,100.,320.):
        assert c.radial_residual(r)<mp.mpf('1e-8')

def test_unsupported_background_is_not_silently_substituted():
    with pytest.raises(ValueError):LiLeaverCloud(a='.877153')
