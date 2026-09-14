"""Local homogeneous gauge jump columns for the 2023 matching construction.

No global radial solution or fitted paper flux enters these columns. This is
only the homogeneous spin-1/kappa block, not the sourced matching RHS.
"""
import numpy as np
from lorenz_jet import Jet
from lorenz_ghp import KerrGHP
from environment_angular_diagnostic import DenseRealHarmonic
from pybhpt.swsh import Yslm


def radial_jet(g,spin,lam,y0,y1):
    n=g.order;y=Jet(y0,n);y.c[1,0]=y1
    K=(g.r*g.r+g.a*g.a)*g.omega-g.a*g.m
    V=(K*K-2j*spin*(g.r-1)*K)/g.delta+4j*spin*g.omega*g.r-lam
    for j in range(n-1):
        rhs=-((spin+1)*2*(g.r-1)*y.derivative(0)+V*y)/g.delta
        y.c[j+2,0]=rhs.c[j,0]/((j+1)*(j+2))
    return y


def angular_jet(g,spin,ell):
    h=DenseRealHarmonic(spin,ell,g.m,g.a*g.omega)
    t=float(g.theta.value.real)
    # separated_jet returns a Kinnersley field S/zeta^(|s|-s), not S.
    # Recur the angular equation directly to avoid importing that tetrad factor.
    out=Jet(h(t),g.order);out.c[0,1]=h(t,deriv=1)
    A=h.eigenvalue-g.a*g.a*g.omega*g.omega+2*g.a*g.m*g.omega
    U=g.a*g.a*g.omega*g.omega*g.theta.cos()**2-2*g.a*g.omega*spin*g.theta.cos()+spin+A-(g.m+spin*g.theta.cos())**2/g.theta.sin()**2
    for j in range(g.order-1):
        rhs=-g.theta.cos()/g.theta.sin()*out.derivative(1)-U*out
        out.c[0,j+2]=rhs.c[0,j]/((j+1)*(j+2))
    return out,h.eigenvalue


def operators(g):
    K=(g.r*g.r+g.a*g.a)*g.omega-g.a*g.m
    Q=g.m/g.theta.sin()-g.a*g.omega*g.theta.sin()
    def D(f,n=0):return f.derivative(0)-1j*K/g.delta*f+n*2*(g.r-1)/g.delta*f
    def Ddag(f,n=0):return f.derivative(0)+1j*K/g.delta*f+n*2*(g.r-1)/g.delta*f
    def L(f,n=0):return f.derivative(1)+Q*f+n*g.theta.cos()/g.theta.sin()*f
    def Ldag(f,n=0):return f.derivative(1)-Q*f+n*g.theta.cos()/g.theta.sin()*f
    return D,Ddag,L,Ldag



def scalar_components(g,h,kappa):
    """Three 2023 trace-sector components, with Box(kappa)=h/2.

    Returned components are h_lp_lp, h_mp_mp, and rho*h_lp_mp.
    Here kappa is the 2023 potential, not the 2024 kappa variable.
    """
    if g.omega==0:raise ValueError('Nonstatic scalar reconstruction only')
    D,_,_,Ldag=operators(g);rho=g.zeta.conjugate()
    return [-(D(g.r*D(h))/(1j*g.omega)+4*D(D(kappa))),
            -Ldag(-g.a*g.theta.cos()/g.omega*Ldag(h)+4*Ldag(kappa),-1),
            -(g.sigma/(2j*g.omega)*D(Ldag(h))
              +g.a/g.omega*(g.r*g.theta.sin()*D(h)+g.theta.cos()*Ldag(h))
              +4*(rho*D(Ldag(kappa))-Ldag(kappa)+1j*g.a*g.theta.sin()*D(kappa)))]


def gauge_basis(g,ell,kind,datum,return_vector=False):
    D,Ddag,L,Ldag=operators(g);rho=g.zeta.conjugate()
    if kind=='spin1':
        Sm,lam=angular_jet(g,-1,ell);Sp,_=angular_jet(g,1,ell)
        Pm=radial_jet(g,-1,lam,1. if datum==0 else 0.,1. if datum==1 else 0.)
        B=np.sqrt(lam*lam+4*g.a*g.m*g.omega-4*g.a*g.a*g.omega*g.omega)
        Pp=g.delta*D(D(Pm))/B;p=(-1)**(ell+g.m)
        calP=D(Pm)+p*Ddag(Pp);calS=Ldag(Sm,1)-p*L(Sp,1)
        if return_vector:
            z=Jet(0.,g.order);one=Jet(1.,g.order)
            lp=[(g.r*g.r+g.a*g.a)/g.delta,one,z,g.a/g.delta]
            lm=[-lp[0],one,z,-lp[3]]
            mp=[1j*g.a*g.theta.sin(),z,one,1j/g.theta.sin()]
            mm=[-mp[0],z,one,-mp[3]]
            return [-(Pm*calS*lp[i]+p*Pp*calS*lm[i]-calP*Sm*mp[i]+p*calP*Sp*mm[i])/(2*g.sigma) for i in range(4)]
        return [p*2/g.delta*D(Pp,-1)*calS,
                p*2*calP*Ldag(Sp,-1),
                p*((rho*D(calP)-2*calP)*Sp+Pp/g.delta*(rho*Ldag(calS)-2*rho.derivative(1)*calS))]
    if kind=='kappa':
        S,lam=angular_jet(g,0,ell)
        k=radial_jet(g,0,lam,1. if datum==0 else 0.,1. if datum==1 else 0.)*S
        if return_vector:return [2*sum(g.inv[i][j]*g.partial(k,j) for j in range(4)) for i in range(4)]
        return [-4*D(D(k)),-4*Ldag(Ldag(k),-1),
                -4*(rho*D(Ldag(k))-Ldag(k)+1j*g.a*g.theta.sin()*D(k))]
    raise ValueError('Unknown homogeneous gauge sector')


def projected_matrix(r0,a,m,ellmax,testmax,quadrature=24):
    if m==0 or ellmax<abs(m) or testmax<abs(m):
        raise ValueError('Nonstatic modes with ellmax,testmax>=|m| required')
    if abs(a)>=1 or r0<=1+np.sqrt(1-a*a):
        raise ValueError('Require subextremal spin and radius outside the horizon')
    x,w=np.polynomial.legendre.leggauss(quadrature);theta=np.arccos(x)
    degrees=list(range(abs(m),testmax+1));labels=[];columns=[]
    spins=(0,2,1)
    angular=np.array([[Yslm(s,j,m,theta) if j>=max(abs(m),abs(s)) else np.zeros_like(theta) for s in spins] for j in degrees])
    for ell in range(abs(m),ellmax+1):
      for kind in ('spin1','kappa'):
       for datum in (0,1):
        samples=[]
        for t in theta:
            g=KerrGHP(r0,t,a,omega=m/(r0**1.5+a),m=m,order=6)
            h=gauge_basis(g,ell,kind,datum)
            samples.append([[v.value for v in h],[v.derivative(0).value for v in h]])
        column=2*np.pi*np.einsum('jfk,kdf,k->djf',angular,np.array(samples),w)
        labels.append((ell,kind,datum));columns.append(column)
    return np.moveaxis(np.array(columns),0,-1),degrees,labels
