"""Static trace-generated Lorenz vacuum piece from arXiv:2306.16459 Sec. IV.

Particular kappa/B solutions have zero data at r0 on each side. Their
homogeneous freedom and the spin-2/completion pieces still require matching.
"""
from functools import lru_cache
import numpy as np
from scipy.integrate import solve_ivp
from environment_source import angular_mode
from lorenz_static_trace import static_trace_radial
from lorenz_ghp import KerrGHP
from lorenz_jet import Jet
from lorenz_tensor import vector_covariant_derivative


@lru_cache(maxsize=32)
def angular_couplings(ell):
    x,w=np.polynomial.legendre.leggauss(32)
    theta=np.arccos(x)
    S,dS,_=angular_mode(theta,ell,0,0.)
    result=[]
    for j in range(max(0,ell-2),ell+3,2):
        T,dT,_=angular_mode(theta,j,0,0.)
        c=2*np.pi*np.dot(w,T*x*x*S)
        d=2*np.pi*np.dot(w,dT*x*x*dS)/(j*(j+1)) if j else 0.
        result.append((j,float(c),float(d)))
    return tuple(result)


@lru_cache(maxsize=128)
def radial_particular(r0,a,ell,side,endpoint):
    modes=angular_couplings(ell)
    lam=ell*(ell+1)
    def rhs(r,y):
        (R,Rp,_),_=static_trace_radial(r,r0,a,ell,side=side if r==r0 else None)
        delta=r*r-2*r+a*a
        values=[]
        for index,(j,c,d) in enumerate(modes):
            k,kp,b,bp=y[4*index:4*index+4]
            sk=R*((r*r if j==ell else 0.)+a*a*c)/2
            sb=delta*Rp*((r*r if j==ell else 0.)+a*a*d)/(2*lam)
            values.extend((kp,(j*(j+1)*k+sk-2*(r-1)*kp)/delta,
                           bp,(j*(j+1)*b+sb)/delta if j else 0.))
        return values
    solution=solve_ivp(rhs,(r0,endpoint),np.zeros(4*len(modes)),method='DOP853',
                       rtol=2e-12,atol=2e-13,dense_output=True)
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution


def _angular_jet(g,ell):
    S,dS,_=angular_mode(float(g.theta.value.real),ell,0,0.)
    out=Jet(S,g.order)
    out.c[0,1]=dS
    for n in range(g.order-1):
        rhs=-g.theta.cos()/g.theta.sin()*out.derivative(1)-ell*(ell+1)*out
        out.c[0,n+2]=rhs.c[0,n]/((n+1)*(n+2))
    return out


def static_trace_metric(r,theta,r0=6.,a=.6,ell=2,order=6):
    if ell<2 or ell%2 or r==r0:
        raise ValueError('Even ell>=2 vacuum points only; monopole completion is separate')
    g=KerrGHP(r,theta,a,omega=0.,m=0,order=order)
    side='inside' if r<r0 else 'outside'
    (R0,Rp,_),_=static_trace_radial(r,r0,a,ell)
    R=Jet(R0,order)
    R.c[1,0]=Rp
    for n in range(order-1):
        rhs=(ell*(ell+1)*R-2*(g.r-1)*R.derivative(0))/g.delta
        R.c[n+2,0]=rhs.c[n,0]/((n+1)*(n+2))
    state=radial_particular(r0,a,ell,side,r).y[:,-1]
    kappa,B=Jet(0.,order),Jet(0.,order)
    for index,(j,c,d) in enumerate(angular_couplings(ell)):
        k0,kp,b0,bp=state[4*index:4*index+4]
        k,b=Jet(k0,order),Jet(b0,order)
        k.c[1,0],b.c[1,0]=kp,bp
        sk=R*((g.r*g.r if j==ell else 0.)+a*a*c)/2
        sb=g.delta*R.derivative(0)*((g.r*g.r if j==ell else 0.)+a*a*d)/(2*ell*(ell+1))
        for n in range(order-1):
            krhs=(j*(j+1)*k+sk-2*(g.r-1)*k.derivative(0))/g.delta
            brhs=(j*(j+1)*b+sb)/g.delta
            k.c[n+2,0]=krhs.c[n,0]/((n+1)*(n+2))
            if j:
                b.c[n+2,0]=brhs.c[n,0]/((n+1)*(n+2))
        S=_angular_jet(g,j)
        kappa=kappa+k*S
        B=B+b*S.derivative(1)
    xi=[Jet(0.,order),kappa.derivative(0)+(B.derivative(1)+g.theta.cos()/g.theta.sin()*B)/g.delta,
        kappa.derivative(1)-B.derivative(0),Jet(0.,order)]
    derivative=vector_covariant_derivative(g,xi)
    h=[[derivative[i][j]+derivative[j][i] for j in range(4)] for i in range(4)]
    return g,h,R*_angular_jet(g,ell)
