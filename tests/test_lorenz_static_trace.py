import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lorenz_static_trace import static_trace_radial


def test_static_trace_equation_and_exact_source_jump():
    for a in (0.,.6,.9):
        for ell in (0,2,4):
            left,J=static_trace_radial(6.,6.,a,ell,'inside')
            right,_=static_trace_radial(6.,6.,a,ell,'outside')
            np.testing.assert_allclose(left[0],right[0],atol=1e-13)
            np.testing.assert_allclose((24+a*a)*(right[1]-left[1]),J,rtol=1e-12)
            for r in (1+np.sqrt(1-a*a),4.,8.,20.):
                (R,dR,ddR),_=static_trace_radial(r,6.,a,ell)
                np.testing.assert_allclose((r*r-2*r+a*a)*ddR+2*(r-1)*dR-ell*(ell+1)*R,0,atol=1e-12)


def test_static_monopole_far_boundary_and_odd_parity():
    value,J=static_trace_radial(1e6,6.,.6,0)
    np.testing.assert_allclose(1e6*value[0],-J,rtol=2e-6)
    for r in (2.,5.,8.):
        odd,_=static_trace_radial(r,6.,.6,3)
        np.testing.assert_allclose(odd,0,atol=1e-14)
