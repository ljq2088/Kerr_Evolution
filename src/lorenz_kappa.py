"""Retarded trace-driven kappa via a mass-squared resolvent derivative.

(Box-lambda) h_lambda=16pi T implies Box d_lambda h_lambda|0=h_0.
Thus kappa=-i omega d_lambda h_lambda|0, with differentiated retarded BCs.
At finite angular truncation the source projector also varies with lambda;
the extra contact terms cancel only in the complete angular sum.
"""
from functools import lru_cache
import numpy as np
from environment_radial import RadialGreen
from environment_source import angular_mode,kerr_metric
from lorenz_mode_jet import separated_jet


@lru_cache(maxsize=96)
def scalar_resolvent(r0,a,ell,m,msq,rtol=1e-12,rmax=2000.):
    op=1/(r0**1.5+a)
    omega=m*op
    metric=kerr_metric(r0,np.pi/2,a)
    ut=1/np.sqrt(-metric[0,0]-2*op*metric[0,3]-op*op*metric[3,3])
    radial=RadialGreen(a,0.,omega,ell,m,rtol=rtol,offset=1e-4,rmax=rmax,mass_squared=msq)
    S=angular_mode(np.pi/2,ell,m,a*a*(omega*omega-msq))[0]
    source=-16*np.pi*S/ut
    zi=radial.insol.sol(r0)[0]*source/radial.w0
    zh=radial.upsol.sol(r0)[0]*source/radial.w0
    return radial,zi,zh


def trace_field_jet(g,r0,ell,msq=0.,rtol=1e-12,rmax=2000.):
    r,theta=float(g.r.value.real),float(g.theta.value.real)
    radial,zi,zh=scalar_resolvent(r0,g.a,ell,g.m,msq,rtol,rmax)
    R,dR=(zh*radial.insol.sol(r) if r<r0 else zi*radial.upsol.sol(r))
    S,dS,lam=angular_mode(theta,ell,g.m,g.a*g.a*(g.omega*g.omega-msq))
    eigenvalue=lam+g.a*g.a*g.omega*g.omega-2*g.a*g.m*g.omega
    return separated_jet(g,0,eigenvalue,R,dR,S,dS,mass_squared=msq)


def kappa_jet(g,r0,ell,step=None,rtol=1e-12,rmax=2000.):
    """Fourth-order central resolvent derivative; vary step/tolerance to audit."""
    if abs(g.omega)<1e-12:
        raise ValueError('Static kappa needs separate treatment')
    if step is None:
        step=min(5e-5,.01*g.omega**2)
    if not 0<step<g.omega**2:
        raise ValueError('Resolvent derivative step must stay below the mass threshold omega^2')
    def h(msq):
        return trace_field_jet(g,r0,ell,msq,rtol,rmax)
    coarse=(h(step)-h(-step))/(2*step)
    fine=(h(step/2)-h(-step/2))/step
    return -1j*g.omega*(4*fine-coarse)/3
