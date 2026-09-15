"""High precision real-frequency angular input for local matching diagnostics.

Assemble the cos(theta) matrix analytically and diagonalize with mpmath.
Evaluate normalized spherical harmonics by factorial Wigner sums. No pybhpt
angular eigenvectors or angular evaluations enter. Outputs retain long double.
"""
from functools import lru_cache
import math
import mpmath as mp
import numpy as np


def number(x):
    return mp.mpf(decimal(x))


def decimal(x):
    # Preserve binary input when passing it to mpmath, including long double.
    return np.format_float_scientific(np.longdouble(x),precision=35,unique=False)


def extended(x):
    return np.longdouble(str(x))


@lru_cache(maxsize=64)
def eigensystem(spin,m,c_text,ellmax=30,dps=50):
    low=max(abs(spin),abs(m));size=ellmax-low+1
    with mp.workdps(dps):
        c=mp.mpf(c_text);X=mp.matrix(size+1)
        for i in range(size+1):
            ell=low+i
            X[i,i]=-mp.mpf(m*spin)/(ell*(ell+1)) if ell else 0
            if i<size:
                j=ell+1
                X[i,i+1]=X[i+1,i]=mp.sqrt(mp.mpf((j*j-m*m)*(j*j-spin*spin))/(4*j*j-1))/j
        square=X*X;H=mp.matrix(size)
        for i in range(size):
            for j in range(size):
                H[i,j]=-c*c*square[i,j]+2*spin*c*X[i,j]
            ell=low+i
            H[i,i]+=ell*(ell+1)-spin*(spin+1)+c*c-2*m*c
        eigen,vectors=mp.eigsy(H)
        return tuple(eigen),tuple(tuple(vectors[i,j]*mp.sign(vectors[j,j]) for i in range(size)) for j in range(size))


@lru_cache(maxsize=65536)
def spherical_pair(spin,ell,m,theta_text,dps=50):
    if ell<max(abs(spin),abs(m)):return np.longdouble(0),np.longdouble(0)
    with mp.workdps(dps):
        t=mp.mpf(theta_text);c=mp.cos(t/2);s=mp.sin(t/2)
        if c==0 or s==0:raise ValueError('Use interior angular points')
        pref=(-1)**spin*mp.sqrt(mp.mpf(2*ell+1)/(4*mp.pi)*math.factorial(ell-spin)*math.factorial(ell+spin)*math.factorial(ell-m)*math.factorial(ell+m))
        value=mp.mpf(0);derivative=mp.mpf(0)
        for k in range(max(0,-m-spin),min(ell-spin,ell-m)+1):
            nc=2*ell-spin-m-2*k;ns=m+spin+2*k
            term=(-1)**(m+spin+k)*c**nc*s**ns/(math.factorial(ell-spin-k)*math.factorial(k)*math.factorial(m+spin+k)*math.factorial(ell-m-k))
            value+=term;derivative+=term*(-nc*s/c+ns*c/s)/2
        return extended(pref*value),extended(pref*derivative)


class PreciseRealHarmonic:
    def __init__(self,spin,ell,m,c,ellmax=30):
        if np.imag(c)!=0:raise ValueError('Real frequency required')
        self.s,self.j,self.m=spin,ell,m;self.jmin=max(abs(spin),abs(m))
        if ell>ellmax-8:raise ValueError('Increase angular basis padding')
        eigen,coeff=eigensystem(spin,m,decimal(np.real(c)),ellmax)
        index=ell-self.jmin
        with mp.workdps(50):
            self.eigenvalue=extended(eigen[index]);self.coeffs=np.array([extended(v) for v in coeff[index]])

    def __call__(self,theta,deriv=0):
        if deriv not in (0,1):raise ValueError('Only value and slope provided; recur the angular ODE for higher derivatives')
        def at(t):
            return sum(b*spherical_pair(self.s,self.jmin+j,self.m,decimal(t))[deriv] for j,b in enumerate(self.coeffs) if abs(b)>1e-42)
        return at(theta) if np.ndim(theta)==0 else np.array([at(t) for t in theta],dtype=np.longdouble)


def precise_angular_jet(g,spin,ell):
    from lorenz_jet import Jet
    h=PreciseRealHarmonic(spin,ell,g.m,g.a*g.omega)
    t=g.theta.value.real
    out=Jet(h(t),g.order);out.c[0,1]=h(t,deriv=1)
    A=h.eigenvalue-g.a*g.a*g.omega*g.omega+2*g.a*g.m*g.omega
    U=g.a*g.a*g.omega*g.omega*g.theta.cos()**2-2*g.a*g.omega*spin*g.theta.cos()+spin+A-(g.m+spin*g.theta.cos())**2/g.theta.sin()**2
    for j in range(g.order-1):
        rhs=-g.theta.cos()/g.theta.sin()*out.derivative(1)-U*out
        out.c[0,j+2]=rhs.c[0,j]/((j+1)*(j+2))
    return out,h.eigenvalue


@lru_cache(maxsize=32)
def precise_trace_angular_data(r0,a,m,L):
    c=a*m/(r0**np.longdouble('1.5')+a)
    eigen,vectors=eigensystem(0,m,decimal(c));low=abs(m);N=len(eigen)
    with mp.workdps(50):
        X=mp.matrix(N+1)
        for i in range(N):
            j=low+i+1;X[i,i+1]=X[i+1,i]=mp.sqrt(mp.mpf(j*j-m*m)/(4*j*j-1))
        C=mp.matrix(N+1,L-low+1)
        for j in range(L-low+1):
            for i in range(N):C[i,j]=vectors[j][i]
        gamma=C.T*X*X*C
        lam=np.array([extended(v) for v in eigen[:L-low+1]])
        eq=np.array([PreciseRealHarmonic(0,ell,m,c)(extended(mp.pi/2)) for ell in range(low,L+1)])
        gamma=np.array([[extended(v) for v in row] for row in gamma.tolist()])
    return lam,eq,gamma


def precise_grid(q):
    with mp.workdps(50):
        x,w=mp.gauss_quadrature(q,'legendre')
        return np.array([extended(mp.acos(v)) for v in x]),np.array([extended(v) for v in w]),extended(mp.pi)


def precise_spherical(spin,ell,m,theta):
    at=lambda t:spherical_pair(spin,ell,m,decimal(t))[0]
    return at(theta) if np.ndim(theta)==0 else np.array([at(t) for t in theta])


def install_matching_angular():
    # Only independent matching modules are changed, not production.
    import paper_jump_basis,paper_sourced_matching,paper_full_tetrad,paper_all_component_matching
    for module in (paper_jump_basis,paper_sourced_matching,paper_full_tetrad,paper_all_component_matching):
        module.angular_jet=precise_angular_jet
    paper_jump_basis.DenseRealHarmonic=PreciseRealHarmonic
    paper_sourced_matching.DenseRealHarmonic=PreciseRealHarmonic
    paper_sourced_matching.trace_angular_data=precise_trace_angular_data
    paper_all_component_matching.trace_angular_data=precise_trace_angular_data
