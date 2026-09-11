"""Kerr static Lorenz completion basis, completion table of 2306.16459v3.

M=1. y is a particular solution with zero radial data at reference_radius;
this harmonic freedom differs from a boundary-adapted paper basis. No particle
matching, boundary regularity or charge normalization is asserted by this module.
"""
from functools import lru_cache
import numpy as np
from scipy.integrate import solve_ivp
from lorenz_ghp import KerrGHP
from lorenz_jet import Jet
from lorenz_tensor import vector_covariant_derivative


@lru_cache(maxsize=64)
def _y_radial(a,reference_radius,r):
    rp,rm=1+np.sqrt(1-a*a),1-np.sqrt(1-a*a)
    def rhs(x,y):
        delta=x*x-2*x+a*a
        f=2/(rp-rm)*np.log((x-rp)/(x-rm))
        return [y[1],((x*x+a*a/3)*f-2*(x-1)*y[1])/delta,
                y[3],(6*y[2]+2*a*a*f/3-2*(x-1)*y[3])/delta]
    if r==reference_radius:
        return np.zeros(4)
    solution=solve_ivp(rhs,(reference_radius,r),np.zeros(4),method='DOP853',rtol=2e-12,atol=2e-14)
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution.y[:,-1]


def completion_scalars(g,reference_radius=6.):
    r,a,c=g.r,g.a,g.theta.cos()
    rp,rm=1+np.sqrt(1-a*a),1-np.sqrt(1-a*a)
    f=2/(rp-rm)*((r-rp)/(r-rm)).log()
    state=_y_radial(a,reference_radius,float(r.value.real))
    modes=[]
    for j,(ell,source) in enumerate(((0,(r*r+a*a/3)*f),(2,2*a*a*f/3))):
        radial=Jet(state[2*j],g.order)
        radial.c[1,0]=state[2*j+1]
        for n in range(g.order-1):
            rhs=(ell*(ell+1)*radial+source-2*(r-1)*radial.derivative(0))/g.delta
            radial.c[n+2,0]=rhs.c[n,0]/((n+1)*(n+2))
        modes.append(radial)
    P2=(3*c*c-1)/2
    y=modes[0]+modes[1]*P2
    z=a/3*(r-rp+2*((r-rm)/(rp-rm)).log()
        -(r-rp)*((rm-5*rp)*r+rp*(rp+3*rm))/(4*(1-a*a))*P2)
    return y,z,f


def _gauge(g,xi):
    d=vector_covariant_derivative(g,xi)
    return [[d[i][j]+d[j][i] for j in range(4)] for i in range(4)]


def completion_metric(r,theta,a=.6,mode='A',order=6,reference_radius=6.):
    if not 0<=a<1 or min(r,reference_radius)<=1+np.sqrt(1-a*a):
        raise ValueError('Exterior subextremal Kerr points required')
    g=KerrGHP(r,theta,a,omega=0.,m=0,order=order)
    r,s,c=g.r,g.theta.sin(),g.theta.cos()
    zero=lambda:Jet(0.,order)
    rp=1+np.sqrt(1-a*a)
    if mode=='A':return g,g.g,Jet(4.,order)
    if mode in ('B','C'):
        xi=[zero(),1/g.delta,zero(),zero()]
        if mode=='B':
            xi=[zero(),(r*(r*r+a*a)-2*rp*rp)/g.delta,a*a*s*c,zero()]
        return g,_gauge(g,xi),Jet(6. if mode=='B' else 0.,order)
    y,z,f=completion_scalars(g,reference_radius)
    if mode in ('D','F'):
        # Lie derivative from xi^r and the explicitly t-linear contravariant
        # components, plus the scalar-gradient gauge part. Result is stationary.
        radial=a*a*c*c/g.sigma if mode=='D' else -a*(r*r-a*a)*c*c/g.sigma
        temporal=[1.,0.,0.,0.] if mode=='D' else [2*a,0.,0.,1.]
        h=[]
        for i in range(4):
            row=[]
            for j in range(4):
                value=radial*g.g[i][j].derivative(0)
                value+=g.g[j][1]*g.partial(radial,i)+g.g[i][1]*g.partial(radial,j)
                if i==0:value+=sum(g.g[j][k]*temporal[k] for k in range(4))
                if j==0:value+=sum(g.g[i][k]*temporal[k] for k in range(4))
                row.append(value)
            h.append(row)
        scalar=y if mode=='D' else z
        extra=_gauge(g,[g.partial(scalar,i) for i in range(4)])
        return g,[[h[i][j]+extra[i][j] for j in range(4)] for i in range(4)],(2+2*f if mode=='D' else Jet(4*a,order))
    h=[[zero() for _ in range(4)] for _ in range(4)]
    if mode=='E':
        h[0][0]=2*r/g.sigma
        h[0][3]=h[3][0]=-2*a*r*s*s/g.sigma
        h[1][1]=2*r*g.sigma/g.delta**2
        h[3][3]=2*a*a*r*s**4/g.sigma
        xi=[g.partial(y,i) for i in range(4)]
        expected=-2*f
    elif mode=='G':
        h[0][0]=-4*a*r*c*c/g.sigma**2
        h[0][3]=h[3][0]=-2*r*s*s/g.sigma+4*a*a*r*s*s*c*c/g.sigma**2
        h[1][1]=2*a*(c*c*g.delta-g.sigma)/g.delta**2
        h[2][2]=2*a*c*c
        h[3][3]=(2*a+4*a*r*s*s/g.sigma-4*a**3*r*s*s*c*c/g.sigma**2)*s*s
        xi=[zero(),a*(r*s*s+c*c)/g.delta,a*s*c,zero()]
        expected=zero()
    else:
        raise ValueError('Completion mode must be A..G')
    extra=_gauge(g,xi)
    return g,[[h[i][j]-extra[i][j] for j in range(4)] for i in range(4)],expected
