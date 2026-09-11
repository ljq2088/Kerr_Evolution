import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lorenz_completion import completion_metric
from lorenz_tensor import trace,lorenz_constraint,linearized_einstein


def test_static_completion_basis_vacuum_trace_and_lorenz():
    for a in (0.,.6):
        for r in (3.5,9.):
            for mode in 'ABCDEFG':
                g,h,expected=completion_metric(r,1.1,a=a,mode=mode,order=6)
                np.testing.assert_allclose(trace(g,h).value,expected.value,rtol=2e-10,atol=2e-10,
                                           err_msg=f'{a=}, {r=}, {mode=}: trace')
                constraint=np.array([v.value for v in lorenz_constraint(g,h)])
                einstein=np.array([[v.value for v in row] for row in linearized_einstein(g,h)])
                np.testing.assert_allclose(constraint,0,atol=2e-9,err_msg=f'{a=}, {r=}, {mode=}: gauge')
                np.testing.assert_allclose(einstein,0,atol=2e-9,err_msg=f'{a=}, {r=}, {mode=}: Einstein')
