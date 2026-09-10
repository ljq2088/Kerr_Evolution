import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lorenz_spin1_chiral import chiral_amplitudes
from lorenz_spin1 import spin1_amplitudes
from lorenz_metric import spin1_metric
from lorenz_tensor import trace,lorenz_constraint,linearized_einstein


def test_self_dual_full_current_matches_contact_shortcut():
    for s in (1,-1):
        parts=chiral_amplitudes(6.,spin=s)
        expected=spin1_amplitudes(6.,spin=s)
        np.testing.assert_allclose([se-dkw for se,dkw in parts],expected,rtol=1e-10,atol=1e-11)
        np.testing.assert_allclose([dkw for se,dkw in parts],0,atol=1e-10)


def test_both_chiralities_preserve_vacuum_lorenz_equations():
    for r in (4.5,8.):
        g,h=spin1_metric(r,1.1,6.,full_current=True)
        scale=max(abs(v.value) for row in h for v in row)
        assert abs(trace(g,h).value)/scale<1e-9
        assert max(abs(v.value) for v in lorenz_constraint(g,h))/scale<1e-9
        assert max(abs(v.value) for row in linearized_einstein(g,h) for v in row)/scale<1e-9
