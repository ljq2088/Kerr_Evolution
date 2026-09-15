import numpy as np
import pytest
from lorenz_ghp import KerrGHP
from lorenz_tensor import vector_covariant_derivative
from paper_jump_basis import gauge_basis
from paper_full_tetrad import vector_tetrad,project_metric,to_bl
from paper_precise_angular import PreciseRealHarmonic
from paper_analytic_jumps import curvature_jumps as analytic_jumps
from paper_sourced_matching import curvature_jumps as ode_jumps
from environment_angular_diagnostic import DenseRealHarmonic


@pytest.mark.parametrize('a,r,m,ell',[(0.,4.,1,1),(.8771530275949366,10.,1,2)])
def test_all_ten_vector_components_are_the_covariant_gauge_transform(a,r,m,ell):
    theta=1.1;g=KerrGHP(r,theta,a,omega=m/(20**1.5+a),m=m,order=8)
    data=(.31+.22j,-.11+.53j)
    bases=[gauge_basis(g,ell,'spin1',datum,True) for datum in (0,1)]
    vector=[sum(data[d]*bases[d][i] for d in (0,1)) for i in range(4)]
    covector=[sum(g.g[i][j]*vector[j] for j in range(4)) for i in range(4)]
    gradient=vector_covariant_derivative(g,covector)
    expected=np.array([[-gradient[i][j].value-gradient[j][i].value for j in range(4)] for i in range(4)])
    actual=vector_tetrad(g,ell,*data)
    projected=project_metric(expected,r,theta,a)
    np.testing.assert_allclose([h.value for h in actual],projected,rtol=2e-10,atol=2e-10)
    np.testing.assert_allclose(to_bl(actual,r,theta,a),expected,rtol=2e-10,atol=2e-10)


@pytest.mark.parametrize('spin,ell,m',[(-2,2,1),(2,6,-1),(0,1,1),(-1,18,2)])
def test_precise_angular_input_matches_independent_double_solver(spin,ell,m):
    c=.013;modes=PreciseRealHarmonic(spin,ell,m,c);reference=DenseRealHarmonic(spin,ell,m,c)
    theta=np.array([.2,.7,1.1,1.8,2.7,3.])
    assert abs(modes.eigenvalue-reference.eigenvalue)<2e-12
    for derivative in (0,1):
        np.testing.assert_allclose(modes(theta,deriv=derivative),reference(theta,deriv=derivative),atol=2e-12,rtol=2e-11)


@pytest.mark.parametrize('a,r,m,ell',[(0.,6.,1,2),(.6,6.,2,2),(.6,6.,-2,3),(.6,6.,3,3),(.8771530275949366,20.,1,18)])
def test_analytic_curvature_jump_matches_independent_radial_amplitudes(a,r,m,ell):
    analytic=np.array(analytic_jumps(r,a,m,ell)[1:])
    reference=np.array(ode_jumps(r,a,m,ell)[1:])
    np.testing.assert_allclose(analytic,reference,rtol=5e-11,atol=2e-11)


def test_algebraic_trace_coupling_matches_direct_sphere_projection():
    from paper_precise_angular import precise_trace_angular_data
    r0=20.;a=.8771530275949366;m=1;L=5;c=a*m/(r0**1.5+a)
    eigen,equator,gamma=precise_trace_angular_data(r0,a,m,L)
    x,w=np.polynomial.legendre.leggauss(80);theta=np.arccos(x)
    reference=[DenseRealHarmonic(0,ell,m,c) for ell in range(m,L+1)]
    samples=np.array([h(theta) for h in reference])
    integrated=2*np.pi*np.einsum('ik,jk,k->ij',samples,samples,w*x*x)
    np.testing.assert_allclose(gamma,integrated,rtol=0,atol=2e-14)
    np.testing.assert_allclose(equator,[h(np.pi/2) for h in reference],rtol=0,atol=2e-14)
    np.testing.assert_allclose(eigen,[h.eigenvalue for h in reference],rtol=0,atol=2e-13)
