import sys
from pathlib import Path
import numpy as np
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
