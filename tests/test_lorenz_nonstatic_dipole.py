import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lorenz_metric import nonstatic_metric
from lorenz_tensor import trace,lorenz_constraint,linearized_einstein
from lorenz_kappa import trace_field_jet


def test_nonstatic_dipole_gauge_sectors_and_reality():
    for r in (4.5,8.):
        g,h=nonstatic_metric(r,1.1,6.,ell=1,m=1,order=8)
        _,negative=nonstatic_metric(r,1.1,6.,ell=1,m=-1,order=8)
        np.testing.assert_allclose(trace(g,h).value,trace_field_jet(g,6.,1).value,rtol=1e-8,atol=1e-9)
        np.testing.assert_allclose([v.value for v in lorenz_constraint(g,h)],0,atol=1e-8)
        np.testing.assert_allclose([[v.value for v in row] for row in linearized_einstein(g,h)],0,atol=1e-8)
        np.testing.assert_allclose([[v.value for v in row] for row in h],
                                  [[v.value.conjugate() for v in row] for row in negative],rtol=1e-7,atol=1e-8)
