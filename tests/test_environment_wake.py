import sys
import json
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from environment_wake import EnvironmentalWake


class ConstantRadial:
    def __init__(self,value):self.value=value
    def evaluate(self,r):return np.array([np.full_like(r,self.value,dtype=complex),np.zeros_like(r)])


def test_helical_phase_axes_and_negative_m_are_independent():
    wake=EnvironmentalWake.__new__(EnvironmentalWake)
    wake.parameters=dict(a=0.,alpha=.3);wake.omega_p=.01;wake.omega_c=.296
    wake.modes={(2,2):(.306,ConstantRadial(1+2j)),(2,-2):(.266,ConstantRadial(3-.5j)),
                (2,0):(.286,ConstantRadial(.5j))}
    theta=np.array([0.,.7,np.pi/2,np.pi]);phi=np.array([.2,.3,.4,.5]);t=17.
    now=wake.evaluate(10.,theta,phi,time=t)
    rotated=wake.evaluate(10.,theta,phi-wake.omega_p*t)*np.exp(-1j*(wake.omega_c-wake.omega_p)*t)
    np.testing.assert_allclose(now,rotated,rtol=1e-13,atol=1e-14)
    # At a=0, S_22=S_2,-2=sqrt(15/(32pi))*sin(theta)^2.
    expected=np.sqrt(15/(32*np.pi))*np.sin(theta)**2*((1+2j)*np.exp(2j*phi)+(3-.5j)*np.exp(-2j*phi))
    expected+=.5j*np.sqrt(5/(16*np.pi))*(3*np.cos(theta)**2-1)
    np.testing.assert_allclose(wake.evaluate(10.,theta,phi),expected,atol=1e-14)


def test_empty_mode_range_is_rejected():
    with pytest.raises(ValueError,match='No modes'):
        EnvironmentalWake([])


def test_one_actual_mode_cannot_be_presented_as_full_wake():
    path=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/forced_mode_nr8_nt10_L6_mg1_sl2_inner0.0005_outer320_log_h32.json'
    report=json.loads(path.read_text())
    with pytest.raises(ValueError,match='Missing scalar modes'):
        EnvironmentalWake([report],ellmax=2)
    candidate=EnvironmentalWake([report],ellmax=2,allow_partial=True)
    assert candidate.status=='partial_mode_field'
    assert candidate.missing==[(2,-2),(2,0)]


def test_actual_modes_with_different_boundaries_require_opt_in():
    directory=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    first=json.loads((directory/'forced_mode_nr8_nt18_L18_mg1_sl2_gh0.0001_go4000_inner0.0005_outer320_log_h32.json').read_text())
    second=json.loads((directory/'forced_mode_nr8_nt18_L18_inner0.0005_outer320_log_h32.json').read_text())
    with pytest.raises(ValueError,match='different discretization'):
        EnvironmentalWake([first,second],ellmax=3,allow_partial=True)
    combined=EnvironmentalWake([first,second],ellmax=3,allow_partial=True,
                                allow_mixed_discretization=True)
    assert combined.mixed_discretization
    assert combined.parameters['green_outer_radius'] is None
    assert combined.parameters['metric_ellmax']==18
    assert combined.radial_domain[1]==1000
    assert [p['numerical']['green_outer_radius'] for p in combined.mode_provenance]==[4000,1000]
    singles=[EnvironmentalWake([row],ellmax=3,allow_partial=True) for row in (first,second)]
    r=np.array([3.,19.,20.,21.,100.])
    np.testing.assert_allclose(combined.evaluate(r,.7,1.2,time=13.),
        sum(wake.evaluate(r,.7,1.2,time=13.) for wake in singles),rtol=1e-13,atol=1e-15)
    with pytest.raises(ValueError,match='common radial solution domain'):
        combined.evaluate(1001.,.7,1.2)
    assert combined.status=='partial_mode_field'


@pytest.mark.parametrize('change',[('alpha',.2),('cloud_mass',.1)])
def test_mixed_grid_opt_in_does_not_allow_different_physics(change):
    import copy
    directory=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    first=json.loads((directory/'forced_mode_nr8_nt18_L18_mg1_sl2_gh0.0001_go4000_inner0.0005_outer320_log_h32.json').read_text())
    second=copy.deepcopy(first)
    second['parameters']['scalar_ell']=4
    second['parameters'][change[0]]=change[1]
    with pytest.raises(ValueError,match='different physics'):
        EnvironmentalWake([first,second],ellmax=4,allow_partial=True,
                          allow_mixed_discretization=True)
