"""Spin-1 asymptotic source amplitudes via the adjoint AAB corrector.

Development implementation: independently compare with public data before use
in physical metric reconstruction. Contact terms at r=r0 are not returned.
"""
import numpy as np
from pybhpt.radial import RadialTeukolsky
from pybhpt.swsh import SpinWeightedSpheroidalHarmonic
from lorenz_ghp import KerrGHP
from lorenz_mode_jet import separated_jet
from lorenz_corrector import corrector


def cky_tensor(g):
    l,n,m,b=g.cov
    kap=[[g.zeta*(m[i]*b[j]-m[j]*b[i]-l[i]*n[j]+l[j]*n[i])/2
          for j in range(4)] for i in range(4)]
    return [[kap[i][j]+kap[i][j].conjugate() for j in range(4)] for i in range(4)]


def current_adjoint_contraction(g,vec,u,fup):
    """u^a u^b times the adjoint j_SE operator, excluding 8pi/ut."""
    dv=[[g.partial(vec[i],j)-sum(g.gamma[k][j][i]*vec[k] for k in range(4))
         for j in range(4)] for i in range(4)]
    sym=[[(dv[i][j]+dv[j][i])/2 for j in range(4)] for i in range(4)]
    nu=corrector(g,sym)
    contraction=sum(u[i]*u[j]*nu[i][j] for i in range(4) for j in range(4))
    trace_term=sum(fup[i][j]*dv[i][j] for i in range(4) for j in range(4))+vec[0]
    return contraction-trace_term


def spin1_amplitudes(r0,a=.6,ell=2,m=2,spin=1):
    if spin not in (-1,1) or m==0:
        raise ValueError('Require nonstatic spin +/-1 mode')
    op=1/(r0**1.5+a)
    omega=m*op
    g=KerrGHP(r0,np.pi/2,a,omega=omega,m=m,order=7)
    metric=np.array([[x.value for x in row] for row in g.g]).real
    ut=1/np.sqrt(-metric[0,0]-2*op*metric[0,3]-op*op*metric[3,3])
    u=np.array([ut,0.,0.,op*ut])
    radial=RadialTeukolsky(spin,ell,m,a,omega,np.array([r0]))
    radial.solve()
    angular=SpinWeightedSpheroidalHarmonic(spin,ell,m,a*omega)
    S0,S1=angular(np.pi/2),angular(np.pi/2,deriv=1)
    ri,di=radial.radialsolution('In',0),radial.radialderivative('In',0)
    ru,du=radial.radialsolution('Up',0),radial.radialderivative('Up',0)
    w0=(r0*r0-2*r0+a*a)**(spin+1)*(ri*du-ru*di)
    f=cky_tensor(g)
    fup=[[sum(g.inv[i][k]*g.inv[j][l]*f[k][l] for k in range(4) for l in range(4))
          for j in range(4)] for i in range(4)]
    amplitudes=[]
    for R0,R1 in ((ri,di),(ru,du)):
        # Fourier phase is reversed for the adjoint kernel; its radial/angular
        # data are the original spin-s mode, with the Sturm-Liouville weight.
        g.omega,g.m=omega,m
        kernel=separated_jet(g,spin,radial.eigenvalue,R0,R1,S0,S1)
        kernel=-2*g.delta**spin*g.zeta**(2*(abs(spin)-spin))*kernel
        g.omega,g.m=-omega,-m
        p=-2*spin
        if spin==1:
            ed=g.derivative(kernel,p,0,'eth')[0]+g.sc['tau']*kernel
            th=g.derivative(kernel,p,0,'thorn')[0]+g.sc['rho']*kernel
            vec=[(-g.cov[0][i]*ed+g.cov[2][i]*th)/2 for i in range(4)]
        else:
            ed=g.derivative(kernel,p,0,'ethp')[0]+g.sc['taup']*kernel
            th=g.derivative(kernel,p,0,'thornp')[0]+g.sc['rhop']*kernel
            vec=[(g.cov[1][i]*ed-g.cov[3][i]*th)/2 for i in range(4)]
        amplitude=8*np.pi*current_adjoint_contraction(g,vec,u,fup).value/(ut*w0)
        amplitudes.append(amplitude)
    return tuple(amplitudes)


if __name__=='__main__':
    for s in (-1,1):
        print(s,spin1_amplitudes(6.,spin=s))
