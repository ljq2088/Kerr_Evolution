import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from paper_jump_basis import gauge_basis
from lorenz_ghp import KerrGHP
from lorenz_tensor import vector_covariant_derivative
from report_nonstatic_paper_projection import transform,PAIRS


def test_explicit_paper_components_against_covariant_gauge_vector():
    for a,r,m,ell in ((0.,6.,1,1),(.8771530275949366,20.,1,1),(.6,6.,2,3)):
      for kind in ('spin1','kappa'):
       for datum in (0,1):
        t=1.1;g=KerrGHP(r,t,a,omega=m/(r**1.5+a),m=m,order=6)
        direct=gauge_basis(g,ell,kind,datum)
        xiup=gauge_basis(g,ell,kind,datum,return_vector=True)
        xicov=[sum(g.g[i][j]*xiup[j] for j in range(4)) for i in range(4)]
        derivative=vector_covariant_derivative(g,xicov)
        metric=[-(derivative[i][j]+derivative[j][i]) for i,j in PAIRS]
        T,dT=transform(r,t,a)
        expected=(T@np.array([h.value for h in metric]))[[0,2,4]]
        dexpected=(T@np.array([h.derivative(0).value for h in metric])+dT@np.array([h.value for h in metric]))[[0,2,4]]
        np.testing.assert_allclose([h.value for h in direct],expected,rtol=2e-10,atol=2e-10)
        np.testing.assert_allclose([h.derivative(0).value for h in direct],dexpected,rtol=2e-10,atol=2e-10)


def test_scalar_convention_map_with_actual_retarded_fields():
    from report_paper_scalar_mapping import check_case
    for case in ((0.,6.,1,1,4.),(.8771530275949366,20.,1,1,10.),(.6,6.,2,2,9.)):
        result=check_case(*case)
        assert result['metric_relative_error']<1e-11
        assert result['derivative_relative_error']<1e-11
        assert result['kappa_equation_relative_residual']<1e-11


def test_angular_jet_is_harmonic_not_kinnersley_field():
    from paper_jump_basis import angular_jet,operators
    from environment_angular_diagnostic import DenseRealHarmonic
    for s in (-2,-1,0,1,2):
        g=KerrGHP(20.,1.1,.8771530275949366,omega=.03,m=1,order=8)
        value,lam=angular_jet(g,s,3)
        reference=DenseRealHarmonic(s,3,1,g.a*g.omega)
        np.testing.assert_allclose(value.value,reference(1.1),rtol=1e-13)
        np.testing.assert_allclose(value.derivative(1).value,reference(1.1,deriv=1),rtol=1e-13)
        assert value.derivative(0).value==0
    D,Dd,L,Ld=operators(g)
    Sm,lam=angular_jet(g,-2,3);Sp,_=angular_jet(g,2,3);aw=g.a*g.omega
    A=np.sqrt(lam**2*(lam+2)**2+8*aw*lam*((g.m-aw)*(5*lam+6)+12*aw)+144*aw**2*(g.m-aw)**2)
    np.testing.assert_allclose(Ld(Ld(Ld(Ld(Sm,2),1)),-1).value,A*Sp.value,rtol=1e-11)


def test_sourced_schwarzschild_jumps_against_published_table():
    from paper_sourced_matching import solve_jumps
    solved=solve_jumps(6.,0.,2,2,10)
    expected=np.array([-269.03958+109.83495j,67.259895-36.611649j,-30.891079,32.178207])
    np.testing.assert_allclose(solved['solution'],expected,rtol=5e-8,atol=1e-6)
    assert solved['unused_residual']<1e-7


def test_direct_paper_kappa_boundary_problem_against_resolvent():
    from report_paper_kappa_radial import compare
    for case in ((0.,6.,1,1,2000.,[3.,5.9,6.1,20.]),
                 (.8771530275949366,20.,1,1,2000.,[1.481,10.,19.9,20.1,40.])):
        row=compare(*case)
        assert row['maximum_relative_kappa_difference']<2e-8
        assert max(r['trace_relative_max'] for r in row['rows'])<2e-8
