import numpy as np
import pytest
from environment_trace_variation import TraceMassVariation
from lorenz_kappa import scalar_resolvent,trace_field_jet,kappa_jet
from environment_source import angular_mode
from lorenz_ghp import KerrGHP


@pytest.mark.parametrize("r0,a,ell,m",[(10,.877,1,1),(20,.6,3,1),(10,.6,2,-2)])
def test_forced_trace_tangent_and_particle_jump(r0,a,ell,m):
    t=TraceMassVariation(r0,a,ell,m,rmax=300)
    r=np.array([3.,r0-.01,r0+.01,50.,150.]);theta=1.1
    step=.001*t.omega**2
    def original(nu):
        g,zi,zh=scalar_resolvent(r0,a,ell,m,nu,rmax=300)
        y=np.where(r<r0,zh*g.insol.sol(r),zi*g.upsol.sol(r))
        s,st,_=angular_mode(theta,ell,m,a*a*(t.omega**2-nu))
        return np.array([y[0]*s,y[1]*s,y[0]*st])
    vals=[original(x*step) for x in (-1,-.5,0,.5,1)]
    fd=(4*(vals[3]-vals[1])/step-(vals[4]-vals[0])/(2*step))/3
    y=t.field_state(r,theta)
    assert np.max(abs(y[:3]-vals[2]))/np.max(abs(vals[2]))<2e-10
    assert np.max(abs(y[3:]-fd))/np.max(abs(fd))<5e-8
    u=t.radial.insol.sol(r0);v=t.radial.upsol.sol(r0)
    jump=t.zi*v[:2]-t.zh*u[:2]
    djump=t.zi_mass*v[:2]+t.zi*v[2:]-t.zh_mass*u[:2]-t.zh*u[2:]
    delta=r0*r0-2*r0+a*a
    np.testing.assert_allclose(jump,[0,t.source/delta],rtol=1e-9,atol=1e-10)
    np.testing.assert_allclose(djump,[0,t.source_mass/delta],rtol=1e-7,atol=1e-10)


def test_eighth_order_jets_and_covariant_kappa_equation():
    t=TraceMassVariation(10,.877,3,1,rmax=300)
    for r in (4.5,12.):
        g=KerrGHP(r,1.1,t.a,omega=t.omega,m=t.m,order=8)
        h,dh=t.field_jet(g);k=t.kappa_jet(g)
        reference=trace_field_jet(g,10,3,rmax=300)
        reference_k=kappa_jet(g,10,3,rmax=300)
        assert np.max(abs(h.c-reference.c))/np.max(abs(reference.c))<2e-10
        assert np.max(abs(k.c-reference_k.c))/np.max(abs(reference_k.c))<5e-8
        np.testing.assert_allclose(g.scalar_wave(k).value,-1j*g.omega*h.value,rtol=2e-11,atol=1e-12)
