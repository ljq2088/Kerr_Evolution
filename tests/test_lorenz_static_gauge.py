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
