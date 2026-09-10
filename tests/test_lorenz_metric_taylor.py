import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lorenz_metric import nonstatic_metric


def test_near_source_taylor_continuation_against_direct_evaluation():
    for sign in (-1.,1.):
        r=6.+sign*5e-5
        dr=-sign*4.5e-5
        _,h=nonstatic_metric(r,1.1,6.,order=8)
        _,direct=nonstatic_metric(r+dr,1.1,6.,order=8)
        values=np.array([[v.value for v in row] for row in direct])
        gradients=np.array([[v.derivative(0).value for v in row] for row in direct])
        approximation=np.array([[(v+dr*v.derivative(0)+dr**2*v.derivative(0).derivative(0)/2).value for v in row] for row in h])
        gradient=np.array([[(v.derivative(0)+dr*v.derivative(0).derivative(0)).value for v in row] for row in h])
        np.testing.assert_allclose(approximation,values,rtol=1e-8,atol=2e-9)
        np.testing.assert_allclose(gradient,gradients,rtol=2e-7,atol=2e-8)
