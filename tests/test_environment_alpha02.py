import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from environment_source import ThresholdCloud


def test_alpha02_cloud_and_mass_normalization():
    cloud=ThresholdCloud(alpha=.2)
    assert 0 < cloud.omega < .2
    e,q=cloud.integrals(4801)
    assert abs(e-1) < 1e-7
    assert abs(e-cloud.omega*q) < 1e-7
    radii=np.geomspace(cloud.rmin,cloud.rmax,1000)
    assert np.all(cloud.radial(radii)[0]>0)  # no radial nodes for |211>
    field,hessian,inv=cloud.hessian(30.,1.2)
    np.testing.assert_allclose(np.einsum('ij,ij',inv,hessian),.2**2*field,rtol=1e-10)
