import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lorenz_static_spin2 import static_spin2_metric,sourced_static_spin2
from lorenz_tensor import trace,lorenz_constraint,linearized_einstein,extreme_weyl


def test_static_hertz_lorenz_vacuum_and_nontrivial_curvature():
    for a,ell,r,branch in ((0.,2,4.5,'P'),(.6,2,8.,'Q'),(.8771530275949366,3,4.5,'P'),(.6,3,8.,'Q')):
        g,h=static_spin2_metric(r,1.1,a=a,ell=ell,branch=branch,order=8)
        scale=max(abs(v.value) for row in h for v in row)
        normalized=[[v/scale for v in row] for row in h]
        assert abs(trace(g,normalized).value)<1e-10
        assert max(abs(v.value) for v in lorenz_constraint(g,normalized))<1e-10
        assert max(abs(v.value) for row in linearized_einstein(g,normalized) for v in row)<1e-10
        _,irg=static_spin2_metric(r,1.1,a=a,ell=ell,branch=branch,order=8,radiation_gauge=True)
        # The complex Hertz reconstruction selects a chirality. The physical
        # static metric includes its conjugate before comparing Weyl scalars.
        reference=np.array(extreme_weyl(g,[[v+v.conjugate() for v in row] for row in irg]))
        physical=np.array(extreme_weyl(g,[[v+v.conjugate() for v in row] for row in h]))
        assert np.max(abs(reference))>1e-10
        np.testing.assert_allclose(physical,reference,rtol=2e-8,atol=1e-11)


def test_static_particle_weyl_inversion_both_parities_and_sides():
    from lorenz_metric import homogeneous_field_jet
    from lorenz_weyl import weyl_amplitudes
    for ell in (2,3):
        for r in (4.5,8.):
            g,h=sourced_static_spin2(r,1.1,ell=ell)
            bc,index=('In',1) if r<6 else ('Up',0)
            amplitudes=weyl_amplitudes(6.,.6,ell,0)
            expected=[homogeneous_field_jet(g,s,ell,amplitudes[s][index],bc).value for s in (2,-2)]
            np.testing.assert_allclose(extreme_weyl(g,h),expected,rtol=2e-9,atol=1e-12)
