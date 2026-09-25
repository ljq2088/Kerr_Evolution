import importlib.util
from pathlib import Path
import numpy as np
import pytest
from pybhpt.swsh import Yslm
from li_separated_source import LiSeparatedSource,SPINS
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('independent_coordinate_fixture',ROOT/'docs/root_cause_followup_20260917/audit_li_appendix_unit_basis.py')
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)

def cloud_angular(t):
 return np.array([np.sin(t)+.17*np.cos(t)**2,
     np.cos(t)-.34*np.sin(t)*np.cos(t),-np.sin(t)-.34*(np.cos(t)**2-np.sin(t)**2)])
def target(t):return np.ones_like(t)/np.sqrt(4*np.pi)

@pytest.mark.parametrize('r',[1.6,3.,20.])
def test_all_ten_components_against_coordinate_hessian(r):
 a=.88;w=.2963;mc=1;mg=0;jmax=3;nq=64
 projection=LiSeparatedSource(a=a,omega_c=w,m_c=mc,metric_m=mg,
     spherical_lmax=jmax,cloud_angular=cloud_angular,target_angular=target,quadrature=nq,pmax=32)
 x,weights=np.polynomial.legendre.leggauss(nq);theta=np.arccos(x)
 for component in range(10):
  j=max(1,abs(SPINS[component]));metric=np.zeros((jmax+1,10),complex)
  metric[j,component]=.31+.73j
  expected=0j;norm=0.
  for t,ww in zip(theta,weights):
   f=audit.fixture(r,t,a,w,mc);weighted=np.zeros(10,complex)
   weighted[component]=metric[j,component]*Yslm(SPINS[component],j,mg,t)
   direct,_=audit.source_targets(weighted,r,t,a,f)
   expected+=2*np.pi*ww*target(t)*direct;norm+=abs(2*np.pi*ww*target(t)*direct)
  state=audit.fixture(r,1.1,a,w,mc)[:3]
  actual=projection.project(r,state,metric)
  assert abs(actual-expected)/max(norm,1e-20)<2e-12,(r,component,actual,expected)

def test_printed_missing_term_changes_mixed_component():
 kw=dict(a=.88,omega_c=.2963,m_c=1,metric_m=1,spherical_lmax=2,
         cloud_angular=cloud_angular,target_angular=target,pmax=32)
 good=LiSeparatedSource(**kw);printed=LiSeparatedSource(**kw,printed_only=True)
 metric=np.zeros((3,10),complex);metric[1,8]=1
 R=audit.fixture(3,1.1,.88,.2963,1)[:3]
 assert abs(good.project(3,R,metric)-printed.project(3,R,metric))>1e-5
