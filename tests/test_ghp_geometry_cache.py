import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lorenz_ghp import KerrGHP,_geometry
from lorenz_jet import Jet


def test_cached_geometry_is_independent_of_fourier_mode_and_mutations():
    _geometry.cache_clear()
    first=KerrGHP(20.,1.1,.877,omega=.02,m=2,order=6)
    expected=first.g[0][0].c.copy()
    first.g[0][0].c[0,0]=999
    first.r.c[0,0]=888
    second=KerrGHP(20.,1.1,.877,omega=-.03,m=-3,order=6)
    assert _geometry.cache_info().hits==1
    np.testing.assert_array_equal(second.g[0][0].c,expected)
    assert second.r.value==20
    assert second.g[0][3] is second.g[3][0]
    field=second.r*second.theta.sin()
    np.testing.assert_array_equal(second.partial(field,0).c,(.03j*field).c)
    np.testing.assert_array_equal(second.partial(field,3).c,(-3j*field).c)


def test_cached_geometry_matches_uncached_coefficients_and_separates_precision(monkeypatch):
    _geometry.cache_clear()
    cached=KerrGHP(8.,.8,.6,order=4)
    direct=KerrGHP.__new__(KerrGHP)
    direct._initialize_geometry(8.,.8,.6,4)
    for i in range(4):
        for j in range(4):
            np.testing.assert_array_equal(cached.g[i][j].c,direct.g[i][j].c)
            for k in range(4):
                np.testing.assert_array_equal(cached.gamma[k][i][j].c,direct.gamma[k][i][j].c)
    monkeypatch.setattr(Jet,'coefficient_dtype',np.clongdouble)
    extended=KerrGHP(8.,.8,.6,order=4)
    assert extended.r.c.dtype==np.dtype(np.clongdouble)
    assert _geometry.cache_info().misses==2
