"""Independent spectral checks of the background and response radial input."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from paper_leaver_cloud import LeaverThresholdCloud
from paper_leaver_radial import LeaverIngoing
from environment_source import ThresholdCloud
from environment_radial import RadialGreen


def test_threshold_series_frequency_profile_and_independent_mass():
    for alpha in (.2,.3):
        leaver=LeaverThresholdCloud(alpha,terms=400)
        old=ThresholdCloud(alpha=alpha);energy,charge=leaver.integrals(128)
        assert abs(float(leaver.a)-old.a)<5e-12
        assert abs(energy/(float(leaver.omega)*charge)-1)<1e-9
        for r in (old.rp+.001,3.,20.,100.,320.):
            np.testing.assert_allclose(leaver.radial(r)/np.sqrt(energy),old.radial(r)[:,0],rtol=1e-7,atol=1e-13)


def test_independent_In_series_convergence_and_horizon_normalization():
    a=.8771530275949366;mu=.3;omega=.28522148715442686
    green=RadialGreen(a,mu,omega,0,0,rmax=1000.,offset=1e-4,rtol=1e-11)
    poor=LeaverIngoing(a,mu,omega,0,0,terms=200)
    fine=LeaverIngoing(a,mu,omega,0,0,terms=3200)
    assert max(abs(poor.state(40.)/green.insol.sol(40.)-1))>.1
    for r in (green.rp+.001,3.,20.,40.):
        value=fine.state(r)
        np.testing.assert_allclose(value,green.insol.sol(r),rtol=1e-8,atol=1e-12)
        charge=-2*(r-green.rp)*(r-(2-green.rp))*np.imag(value[0].conjugate()*value[1])
        np.testing.assert_allclose(charge,4*green.rp*omega,rtol=1e-9)
