"""Development calculation of the compact-source auxiliary Lorenz scalar.

The extended trace-driven kappa term must be treated separately. Do not regard
this routine as a full reconstruction until its normalization is audited.
"""
import numpy as np
from pybhpt.radial import RadialTeukolsky
from pybhpt.swsh import SpinWeightedSpheroidalHarmonic
from lorenz_ghp import KerrGHP
from lorenz_mode_jet import separated_jet
from lorenz_corrector import Weighted
from lorenz_spin1 import cky_tensor,current_adjoint_contraction
from functools import lru_cache


def sdag_tensor(g,f,which):
    """Gravitational S_0^dag or S_4^dag, Appendix B.2, covariant BL output."""
    if which==0:
        p,l,m= -4,0,2
        dk,ek,rk,tk='thorn','eth','rho','tau'
        tb=g.sc['taup'].conjugate()
    else:
        p,l,m=4,1,3
        dk,ek,rk,tk='thornp','ethp','rhop','taup'
        tb=g.sc['tau'].conjugate()
    F=Weighted(f,p,0)
    shift=1 if which==0 else -1
    rho=Weighted(g.sc[rk],shift,shift)
    tau=Weighted(g.sc[tk],shift,-shift)
    bar_rho=rho.conjugate()
    bar_tau_other=Weighted(tb,shift,-shift)
    def D(x,k):
        v,p,q=g.derivative(x.f,x.p,x.q,k)
        return Weighted(v,p,q)
    d=lambda x:D(x,dk)
    e=lambda x:D(x,ek)
    ei=e(F)+3*tau*F
    di=d(F)+3*rho*F
    A=-(e(ei)-tau*ei)/2
    B=-(d(di)-rho*di)/2
    C=(d(ei)+(-rho+bar_rho)*ei+e(di)+(-tau+bar_tau_other)*di)/2
    return [[g.cov[l][a]*g.cov[l][b]*A.f+g.cov[m][a]*g.cov[m][b]*B.f
             +(g.cov[l][a]*g.cov[m][b]+g.cov[m][a]*g.cov[l][b])*C.f/2
             for b in range(4)] for a in range(4)]


def chi_adjoint(g,kernel,which):
    """chi_0^dag or chi_4^dag acting on a scalar adjoint Green function."""
    kinds=('ethp','thornp') if which==0 else ('eth','thorn')
    F=Weighted(kernel)
    def adj(x,k):
        v,p,q=g.derivative(g.sigma*x.f,x.p,x.q,k)
        return Weighted(-v/g.sigma,p,q)
    for _ in range(2):
        F=adj(F,kinds[0])
    F=Weighted(g.zeta.conjugate()**2*F.f,F.p,F.q)
    for _ in range(2):
        F=adj(F,kinds[1])
    return F.f


def direct_chi(g,kernel,which,p=0,q=0):
    kinds=('thornp','ethp') if which==0 else ('thorn','eth')
    F=Weighted(kernel,p,q)
    for _ in range(2):
        v,p,q=g.derivative(F.f,F.p,F.q,kinds[0])
        F=Weighted(v,p,q)
    F=Weighted(g.zeta.conjugate()**2*F.f,F.p,F.q)
    for _ in range(2):
        v,p,q=g.derivative(F.f,F.p,F.q,kinds[1])
        F=Weighted(v,p,q)
    return F.f


def dkw_current_adjoint(g,vec):
    """Adjoint of the explicitly displayed j_DKWSE, without 4/(9 Lt)."""
    # vec is covariant; contractions with tetrads give V^a e_a.
    def derivative(x,k,adjoint=False):
        arg=g.sigma*x.f if adjoint else x.f
        v,p,q=g.derivative(arg,x.p,x.q,k)
        return Weighted(-v/g.sigma if adjoint else v,p,q)
    outputs=[]
    for which,(A,B,dk,ek,pA,qA,pB,qB) in (
            (0,(1,3,'thornp','ethp',-1,-1,-1,1)),
            (4,(0,2,'thorn','eth',1,1,1,-1))):
        wa=Weighted(sum(g.tetrad[A][i]*vec[i] for i in range(4))/g.zeta,pA,qA)
        wb=Weighted(sum(g.tetrad[B][i]*vec[i] for i in range(4))/g.zeta,pB,qB)
        F=derivative(wa,ek,True)-derivative(wb,dk,True)
        F=Weighted(g.zeta*F.f,F.p,F.q)
        F=derivative(F,dk)
        F=Weighted(g.zeta.conjugate()*F.f,F.p,F.q)
        F=derivative(F,ek)
        outputs.append(sdag_tensor(g,g.zeta**4*F.f,which))
    return [[outputs[0][i][j]+outputs[1][i][j] for j in range(4)] for i in range(4)]


@lru_cache(maxsize=256)
def chi_amplitudes(r0,a=.6,ell=2,m=2,return_parts=False):
    if m==0:
        raise ValueError('Nonstatic modes only')
    op=1/(r0**1.5+a)
    omega=m*op
    g=KerrGHP(r0,np.pi/2,a,omega=omega,m=m,order=8)
    metric=np.array([[v.value for v in row] for row in g.g]).real
    ut=1/np.sqrt(-metric[0,0]-2*op*metric[0,3]-op*op*metric[3,3])
    u=np.array([ut,0.,0.,op*ut])
    radial=RadialTeukolsky(0,ell,m,a,omega,np.array([r0]))
    radial.solve()
    angular=SpinWeightedSpheroidalHarmonic(0,ell,m,a*omega)
    S0,S1=angular(np.pi/2),angular(np.pi/2,deriv=1)
    ri,di=radial.radialsolution('In',0),radial.radialderivative('In',0)
    ru,du=radial.radialsolution('Up',0),radial.radialderivative('Up',0)
    w0=(r0*r0-2*r0+a*a)*(ri*du-ru*di)
    f=cky_tensor(g)
    fup=[[sum(g.inv[i][k]*g.inv[j][l]*f[k][l] for k in range(4) for l in range(4))
          for j in range(4)] for i in range(4)]
    amplitudes=[]
    parts=[]
    for R0,R1 in ((ri,di),(ru,du)):
        g.omega,g.m=omega,m
        kernel=separated_jet(g,0,radial.eigenvalue,R0,R1,S0,S1)
        g.omega,g.m=-omega,-m
        vec=[sum(g.g[i][b]*fup[c][b]*g.partial(kernel,c) for b in range(4) for c in range(4))
             for i in range(4)]
        first=-1j/omega*current_adjoint_contraction(g,vec,u,fup)
        p4=sdag_tensor(g,g.zeta**4*direct_chi(g,kernel,4),4)
        p0=sdag_tensor(g,g.zeta**4*direct_chi(g,kernel,0),0)
        current=dkw_current_adjoint(g,vec)
        third=-sum(u[i]*u[j]*(p4[i][j]-p0[i][j]+4*current[i][j])
                   for i in range(4) for j in range(4))/(9*omega**2)
        amplitudes.append(8*np.pi*(first+third).value/(ut*w0))
        parts.append((8*np.pi*first.value/(ut*w0),8*np.pi*third.value/(ut*w0)))
    return parts if return_parts else tuple(amplitudes)


if __name__=='__main__':
    for r in (4.,6.):
        print(r,chi_amplitudes(r,return_parts=True))
