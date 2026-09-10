import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lorenz_ghp import KerrGHP
from lorenz_tensor import vector_covariant_derivative,linearized_einstein,lorenz_constraint,trace
from lorenz_metric import spin2_metric,spin1_metric,spin0_metric
from lorenz_kappa import trace_field_jet
from lorenz_spin1 import cky_tensor,spin1_amplitudes
from lorenz_tensor import tensor_divergence,extreme_weyl
from lorenz_metric import homogeneous_field_jet
from lorenz_weyl import weyl_amplitudes


def test_linearized_einstein_annihilates_pure_gauge():
    g=KerrGHP(8.,1.1,.6,omega=.13,m=2,order=5)
    xi=[g.r*g.theta.sin(),g.r**2,g.theta.cos(),g.r*g.theta.cos()]
    d=vector_covariant_derivative(g,xi)
    h=[[d[i][j]+d[j][i] for j in range(4)] for i in range(4)]
    result=linearized_einstein(g,h)
    np.testing.assert_allclose([[x.value for x in row] for row in result],0,atol=2e-11)
    np.testing.assert_allclose(extreme_weyl(g,h),0,atol=2e-11)


def test_circularity_inversion_coefficient():
    g=KerrGHP(8.,1.1,.6,omega=.13,m=2,order=5)
    xi=[g.r*g.theta.sin(),g.r**2,g.theta.cos(),g.r*g.theta.cos()]
    d=vector_covariant_derivative(g,xi)
    F=[[d[j][i]-d[i][j] for j in range(4)] for i in range(4)]
    f=cky_tensor(g)
    raised=[[sum(f[i][j]*g.inv[j][k] for j in range(4)) for k in range(4)] for i in range(4)]
    H=[[sum((raised[i][k]*F[k][j]-raised[j][k]*F[k][i])/2 for k in range(4))
        for j in range(4)] for i in range(4)]
    divH=[-x for x in tensor_divergence(g,H)]
    current=tensor_divergence(g,F)
    scalar=xi[0]-sum(g.inv[a][i]*g.inv[b][j]*f[i][j]*F[a][b]/2
                    for a in range(4) for b in range(4) for i in range(4) for j in range(4))
    for a in range(4):
        recovered=2*divH[a]+sum(raised[a][b]*current[b] for b in range(4))+g.partial(scalar,a)
        np.testing.assert_allclose(recovered.value,(-1j*g.omega*xi[a]).value,atol=2e-11)


def test_spin2_vacuum_metric():
    for r in (4.5,8.):
        g,h=spin2_metric(r,1.1,6.)
        scale=max(abs(v.value) for row in h for v in row)
        tr=abs(trace(g,h).value)
        gauge=max(abs(v.value) for v in lorenz_constraint(g,h))
        einstein=max(abs(v.value) for row in linearized_einstein(g,h) for v in row)
        print(r,scale,tr,gauge,einstein)
        assert gauge/scale<1e-8
        assert tr/scale<1e-8
        assert einstein/scale<1e-8
        bc,index=('In',1) if r<6 else ('Up',0)
        amps=weyl_amplitudes(6.,.6,2,2)
        expected=[homogeneous_field_jet(g,s,2,amps[s][index],bc).value for s in (2,-2)]
        np.testing.assert_allclose(extreme_weyl(g,h),expected,rtol=1e-8,atol=1e-12)


def test_spin0_vacuum_metric_and_trace():
    for r in (4.5,8.):
        g,h=spin0_metric(r,1.1,6.)
        scale=max(abs(v.value) for row in h for v in row)
        expected=trace_field_jet(g,6.,2).value
        tr=trace(g,h).value
        gauge=max(abs(v.value) for v in lorenz_constraint(g,h))
        einstein=max(abs(v.value) for row in linearized_einstein(g,h) for v in row)
        print('spin0',r,scale,tr,expected,gauge,einstein)
        np.testing.assert_allclose(tr,expected,rtol=1e-7,atol=1e-9)
        assert gauge/scale<1e-7
        assert einstein/scale<1e-8


def test_spin1_vacuum_metric():
    for r in (4.5,8.):
        g,h=spin1_metric(r,1.1,6.)
        scale=max(abs(v.value) for row in h for v in row)
        tr=abs(trace(g,h).value)
        gauge=max(abs(v.value) for v in lorenz_constraint(g,h))
        einstein=max(abs(v.value) for row in linearized_einstein(g,h) for v in row)
        print('spin1',r,scale,tr,gauge,einstein)
        assert gauge/scale<1e-8
        assert tr/scale<1e-8
        assert einstein/scale<1e-8


def test_spin1_recovers_input_maxwell_scalars():
    for r in (4.5,8.):
        g,xi=spin1_metric(r,1.1,6.,return_vector=True)
        d=vector_covariant_derivative(g,xi)
        F=[[d[j][i]-d[i][j] for j in range(4)] for i in range(4)]
        l,n,mm,mb=g.tetrad
        computed=[sum(F[i][j]*u[i]*v[j] for i in range(4) for j in range(4)).value
                  for u,v in ((l,mm),(mb,n))]
        bc,index=('In',1) if r<6 else ('Up',0)
        expected=[homogeneous_field_jet(g,s,2,spin1_amplitudes(6.,.6,2,2,s)[index],bc).value for s in (1,-1)]
        np.testing.assert_allclose(computed,expected,rtol=1e-9,atol=1e-12)
