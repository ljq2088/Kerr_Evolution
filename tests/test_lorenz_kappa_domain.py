"""Trace/kappa must reject dense-output extrapolation at the public entry point."""
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import lorenz_kappa


def point(r):
    return SimpleNamespace(r=SimpleNamespace(value=complex(r)),
                           theta=SimpleNamespace(value=1.1+0j),
                           a=.6,m=1,omega=.1)


@pytest.fixture
def radial_stub(monkeypatch):
    calls=[]
    radial_values=np.array([1+2j,3-5j])
    def solution(branch):
        def evaluate(r):
            calls.append((branch,r))
            return radial_values.copy()
        return SimpleNamespace(sol=evaluate)
    radial=SimpleNamespace(rmin=2.,rmax=8.,insol=solution('In'),upsol=solution('Up'))
    zi,zh=2-1j,-3+4j
    monkeypatch.setattr(lorenz_kappa,'scalar_resolvent',lambda *a,**kw:(radial,zi,zh))
    monkeypatch.setattr(lorenz_kappa,'angular_mode',lambda *a:(.7,.2,4.))
    monkeypatch.setattr(lorenz_kappa,'separated_jet',
                        lambda g,s,eigenvalue,R,dR,S,dS,**kw:np.array([R,dR]))
    return calls,radial_values,zi,zh


@pytest.mark.parametrize('r',[np.nextafter(2.,-np.inf),np.nextafter(8.,np.inf),
                             np.nan,np.inf,-np.inf])
@pytest.mark.parametrize('evaluate',[lorenz_kappa.trace_field_jet,lorenz_kappa.kappa_jet])
def test_trace_and_kappa_reject_invalid_radius_before_dense_output(radial_stub,r,evaluate):
    calls,*_=radial_stub
    with pytest.raises(ValueError) as error:
        evaluate(point(r),5.,1)
    message=str(error.value)
    assert f'requested r={float(r)!r}' in message
    assert 'solved radial range [2.0, 8.0]' in message
    assert not calls


@pytest.mark.parametrize('r,branch',[(2.,'In'),(3.,'In'),(5.,'Up'),(7.,'Up'),(8.,'Up')])
def test_trace_includes_endpoints_and_preserves_valid_branch_values(radial_stub,r,branch):
    calls,radial_values,zi,zh=radial_stub
    actual=lorenz_kappa.trace_field_jet(point(r),5.,1)
    np.testing.assert_array_equal(actual,(zh if branch=='In' else zi)*radial_values)
    assert calls==[(branch,r)]
