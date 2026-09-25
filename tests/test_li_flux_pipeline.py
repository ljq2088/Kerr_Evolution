import numpy as np
import pytest
from li_normalized_cloud import LiNormalizedCloud
from li_order4_green import LiOrder4Green
from environment_radial import RadialGreen
from environment_kerr_quasibound_cloud import GeneralKerrQuasiboundCloud

@pytest.fixture(scope='module')
def cloud():return LiNormalizedCloud()

def test_vector_leaver_evaluation_matches_arbitrary_precision(cloud):
    for r in (cloud.rp+5e-4,2.,20.,320.,1000.):
        value=cloud.radial_state(r)[:,0]/cloud.amplitude
        ref=np.array([complex(cloud.leaver.radial_mp(r,k)) for k in range(3)])
        assert np.max(abs(value/ref-1))<2e-9

def test_bl_unit_mass_and_cutoff_control(cloud):
    assert abs(cloud.bl_integrals(1e-6,2000.)['energy']-1)<2e-13
    assert max(abs(v['energy']-1) for v in cloud.normalization_audit())<1e-6
    ref=GeneralKerrQuasiboundCloud(freeze_growth=True)
    for r in (3.,20.,100.):
        # Independent shooting and a different slice normalization.
        assert abs(abs(cloud.radial(r)[0,0]/ref.radial(r)[0,0])-1)<2e-6

@pytest.mark.parametrize('ell,m',[(0,0),(2,2)])
def test_order4_response_against_independent_coulomb_boundary(ell,m):
    a=.88;mu=.3;omega=.29629353472811227+(m-1)/(20**1.5+a)
    g=LiOrder4Green(a,mu,omega,ell,m)
    ref=RadialGreen(a,mu,omega,ell,m,rmax=g.rmax,offset=1e-4,rtol=1e-11,infinity_method='coulomb')
    r=np.geomspace(2.,320.,121)
    # Compare Green kernels, invariant under either homogeneous rescaling.
    k=g.insol.sol(np.array([3.]))[0,0]*g.upsol.sol(r)[0]/g.w0
    kr=ref.insol.sol(np.array([3.]))[0,0]*ref.upsol.sol(r)[0]/ref.w0
    assert np.linalg.norm(k-kr)/np.linalg.norm(kr)<2e-6
    assert np.max(abs(g.wronskian(r)/g.w0-1))<1e-6
