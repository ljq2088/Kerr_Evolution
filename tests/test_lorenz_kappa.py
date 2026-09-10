import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lorenz_ghp import KerrGHP
from lorenz_kappa import trace_field_jet,kappa_jet


def test_retarded_resolvent_kappa_equation():
    for r in (4.5,8.):
        g=KerrGHP(r,1.1,.6,omega=2/(6**1.5+.6),m=2,order=6)
        h=trace_field_jet(g,6.,2)
        k=kappa_jet(g,6.,2)
        target=-1j*g.omega*h.value
        computed=g.scalar_wave(k).value
        print(r,computed,target,abs((computed-target)/target))
        np.testing.assert_allclose(computed,target,rtol=2e-6,atol=1e-10)


def test_low_frequency_resolvent_stays_on_radiative_branch():
    a=.8771530275949366
    omega=2/(41.8**1.5+a)
    g=KerrGHP(10.,1.1,a,omega=omega,m=2,order=6)
    value=kappa_jet(g,41.8,2).value
    finer=kappa_jet(g,41.8,2,step=.005*omega**2).value
    np.testing.assert_allclose(value,finer,rtol=2e-5,atol=1e-8)
    with pytest.raises(ValueError,match='mass threshold'):
        kappa_jet(g,41.8,2,step=omega**2)
