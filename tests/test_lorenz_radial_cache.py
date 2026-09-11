import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from pybhpt.radial import RadialTeukolsky
from lorenz_ghp import KerrGHP
from lorenz_metric import _homogeneous_radial_data,homogeneous_field_jet


def test_radial_cache_preserves_direct_solver_values_and_boundary_keys():
    _homogeneous_radial_data.cache_clear()
    for spin,ell,m,a,omega,r,bc in ((2,2,2,.6,.13,8.,'Up'),
                                  (2,2,2,.6,.13,8.,'In'),
                                  (-1,6,-2,.877,-.022,30.,'Up')):
        radial=RadialTeukolsky(spin,ell,m,a,omega,np.array([r]));radial.solve(bc=bc)
        expected=(radial.eigenvalue,radial.radialsolution(bc,0),radial.radialderivative(bc,0))
        actual=_homogeneous_radial_data(spin,ell,m,a,omega,r,bc)
        np.testing.assert_array_equal(actual,expected)
        assert all(np.isscalar(value) for value in actual)
    assert _homogeneous_radial_data.cache_info().misses==3


def test_angles_and_amplitudes_reuse_radial_data(monkeypatch):
    import lorenz_metric
    from environment_angular_diagnostic import DenseRealHarmonic
    monkeypatch.setattr(lorenz_metric,'SpinWeightedSpheroidalHarmonic',DenseRealHarmonic)
    _homogeneous_radial_data.cache_clear()
    first=KerrGHP(8.,.7,.6,omega=.13,m=2,order=6)
    second=KerrGHP(8.,1.1,.6,omega=.13,m=2,order=6)
    homogeneous_field_jet(first,2,2,1.,'Up')
    value=homogeneous_field_jet(second,2,2,1.,'Up')
    scaled=homogeneous_field_jet(second,2,2,2.,'Up')
    assert _homogeneous_radial_data.cache_info().misses==1
    assert _homogeneous_radial_data.cache_info().hits==2
    np.testing.assert_array_equal(scaled.c,2*value.c)
