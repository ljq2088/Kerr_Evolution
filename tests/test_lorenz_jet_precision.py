import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lorenz_jet import Jet


def test_extended_jet_retains_small_polynomial_coefficients(monkeypatch):
    if np.finfo(np.longdouble).eps>=np.finfo(float).eps:
        pytest.skip('Extended precision is unavailable')
    monkeypatch.setattr(Jet,'coefficient_dtype',np.clongdouble)
    x=Jet.variable(np.longdouble('0.3'),0,6)
    y=Jet.variable(np.longdouble('0.4'),1,6)
    # Independent polynomial identity, with subtraction that double cannot resolve.
    tiny=np.longdouble('1e-18')
    result=((1+x*y)+tiny*(x+y))-(1+x*y)
    assert abs(result.c[1,0]-tiny)<np.longdouble('2e-20')
    assert abs(result.c[0,1]-tiny)<np.longdouble('2e-20')
    assert result.c.dtype==np.dtype(np.clongdouble)
    identity=(1+x+y)*(1/(1+x+y))
    expected=np.zeros_like(identity.c); expected[0,0]=1
    np.testing.assert_allclose(identity.c,expected,rtol=0,atol=2e-17)
