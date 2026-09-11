import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lorenz_static_gauge import static_trace_metric
from lorenz_tensor import trace,lorenz_constraint,linearized_einstein


def test_static_trace_metric_lorenz_and_vacuum_equations():
    for a,r in ((0.,4.5),(.6,4.5),(.6,8.)):
        g,h,expected=static_trace_metric(r,1.1,a=a)
        np.testing.assert_allclose(trace(g,h).value,expected.value,rtol=1e-9,atol=1e-10)
        np.testing.assert_allclose([v.value for v in lorenz_constraint(g,h)],0,atol=1e-9)
        np.testing.assert_allclose([[v.value for v in row] for row in linearized_einstein(g,h)],0,atol=1e-9)
def test_high_degree_trace_couplings_match_legendre_recurrence():
    from lorenz_static_gauge import angular_couplings, _angular_jet
    from lorenz_ghp import KerrGHP
    from scipy.special import eval_legendre
    import numpy as np
    for ell in (18,20):
        c=lambda n:n/np.sqrt((2*n-1)*(2*n+1))
        expected={ell-2:c(ell)*c(ell-1),ell:c(ell)**2+c(ell+1)**2,
                  ell+2:c(ell+1)*c(ell+2)}
        for j,scalar,gradient in angular_couplings(ell):
            np.testing.assert_allclose(scalar,expected[j],rtol=1e-12,atol=1e-14)
            analytic=((j==ell)+(.5*(ell*(ell+1)+j*(j+1))-3)*expected[j])/(j*(j+1))
            np.testing.assert_allclose(gradient,analytic,rtol=1e-12,atol=1e-14)
        g=KerrGHP(20.,1.1,.8771530275949366,order=4)
        jet=_angular_jet(g,ell+2)
        np.testing.assert_allclose(jet.value,np.sqrt((2*(ell+2)+1)/(4*np.pi))*
                                   eval_legendre(ell+2,np.cos(1.1)),rtol=1e-12)
