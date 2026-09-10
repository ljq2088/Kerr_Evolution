"""Exact m=omega=0 trace Green function for a circular Kerr particle.

This solves only Box h=16 pi T. It is not a static Lorenz metric completion.
M=particle mass=1; spherical harmonics have unit norm over the sphere.
"""
import numpy as np
from scipy.special import eval_legendre,lqmn
from environment_source import angular_mode,kerr_metric


def static_trace_radial(r,r0,a=.6,ell=0,side=None):
    """Return R,R',R'' for h_l0=R(r)Y_l0; select side at the particle."""
    if not 0<=abs(a)<1 or ell<0 or int(ell)!=ell:
        raise ValueError('Require subextremal Kerr and integer ell>=0')
    b=np.sqrt(1-a*a)
    rp=1+b
    if r<rp or r0<=rp or (r==r0 and side not in ('inside','outside')):
        raise ValueError('Exterior points required; choose a side at r=r0')
    op=1/(r0**1.5+a)
    g=kerr_metric(r0,np.pi/2,a)
    ut=1/np.sqrt(-g[0,0]-2*op*g[0,3]-op*op*g[3,3])
    if not np.isfinite(ut):
        raise ValueError('Orbit must be timelike')
    S=angular_mode(np.pi/2,ell,0,0.)[0]
    J=-16*np.pi*S/ut
    x,x0=(r-1)/b,(r0-1)/b
    inside=(r<r0 or (r==r0 and side=='inside'))
    if inside:
        q0=lqmn(0,ell,x0)[0][0,ell]
        R=J/(-b)*eval_legendre(ell,x)*q0
        # Polynomial derivative remains regular at the horizon x=1.
        polynomial=np.polynomial.legendre.Legendre.basis(ell)
        dR=J/(-b)*polynomial.deriv()(x)*q0/b
        ddR=J/(-b)*polynomial.deriv(2)(x)*q0/b**2
    else:
        q,dq=lqmn(0,ell,x)
        coefficient=J/(-b)*eval_legendre(ell,x0)
        R=coefficient*q[0,ell]
        dR=coefficient*dq[0,ell]/b
        ddR=(ell*(ell+1)*R-2*(r-1)*dR)/(r*r-2*r+a*a)
    return np.array([R,dR,ddR]),J
