import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from report_static_matching import scalar_hessian
from lorenz_static_spin2 import sourced_static_spin2
from lorenz_tensor import trace,lorenz_constraint,linearized_einstein,extreme_weyl


def test_free_scalar_hessian_is_traceless_lorenz_vacuum():
    for a,ell,datum in ((0.,2,0),(.6,2,1),(.6,4,0)):
        g,h=scalar_hessian(6.,1.1,a,ell,datum)
        np.testing.assert_allclose(trace(g,h).value,0,atol=2e-11)
        np.testing.assert_allclose([v.value for v in lorenz_constraint(g,h)],0,atol=2e-11)
        np.testing.assert_allclose([[v.value for v in row] for row in linearized_einstein(g,h)],0,atol=2e-11)


def test_circular_isometry_preserves_static_curvature_and_vacuum():
    for a,r,ell in ((0.,4.5,3),(.6,4.5,2),(.6,8.,3)):
        g,raw=sourced_static_spin2(r,1.1,a=a,ell=ell,circular_symmetry=False)
        _,h=sourced_static_spin2(r,1.1,a=a,ell=ell)
        np.testing.assert_allclose(extreme_weyl(g,h),extreme_weyl(g,raw),rtol=2e-10,atol=1e-12)
        np.testing.assert_allclose([v.value for v in lorenz_constraint(g,h)],0,atol=2e-10)
        np.testing.assert_allclose([[v.value for v in row] for row in linearized_einstein(g,h)],0,atol=2e-10)
        for i,j in ((0,1),(0,2),(1,3),(2,3)):
            assert h[i][j].value==0
