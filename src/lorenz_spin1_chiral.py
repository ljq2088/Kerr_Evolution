"""Both Maxwell chiralities from the full complex compact current.

Unlike the self-dual shortcut, retain the explicit DKW source contribution.
This is a development source audit, not yet a production metric interface.
"""
from functools import lru_cache
import numpy as np
from pybhpt.radial import RadialTeukolsky
from pybhpt.swsh import SpinWeightedSpheroidalHarmonic
from lorenz_ghp import KerrGHP
from lorenz_mode_jet import separated_jet
from lorenz_spin1 import cky_tensor,current_adjoint_contraction
from lorenz_chi import dkw_current_adjoint


@lru_cache(maxsize=256)
def chiral_amplitudes(r0,a=.6,ell=2,m=2,spin=1,chirality=1):
    if chirality not in (-1,1) or spin not in (-1,1) or m==0:
        raise ValueError('Nonstatic Maxwell modes and chirality +/-1 required')
    op=1/(r0**1.5+a)
    omega=m*op
    mf,wf=chirality*m,chirality*omega
    g=KerrGHP(r0,np.pi/2,a,omega=wf,m=mf,order=8)
    metric=np.array([[v.value for v in row] for row in g.g]).real
    ut=1/np.sqrt(-metric[0,0]-2*op*metric[0,3]-op*op*metric[3,3])
    u=np.array([ut,0.,0.,op*ut])
    radial=RadialTeukolsky(spin,ell,mf,a,wf,np.array([r0]))
    radial.solve()
    angular=SpinWeightedSpheroidalHarmonic(spin,ell,mf,a*wf)
    ri,di=radial.radialsolution('In',0),radial.radialderivative('In',0)
    ru,du=radial.radialsolution('Up',0),radial.radialderivative('Up',0)
    w0=(r0*r0-2*r0+a*a)**(spin+1)*(ri*du-ru*di)
    if chirality<0:
        w0=w0.conjugate()
    f=cky_tensor(g)
    fup=[[sum(g.inv[i][k]*g.inv[j][l]*f[k][l] for k in range(4) for l in range(4)) for j in range(4)] for i in range(4)]
    results=[]
    for R0,R1 in ((ri,di),(ru,du)):
        g.omega,g.m=wf,mf
        kernel=separated_jet(g,spin,radial.eigenvalue,R0,R1,angular(np.pi/2),angular(np.pi/2,deriv=1))
        kernel=-2*g.delta**spin*g.zeta**(2*(abs(spin)-spin))*kernel
        g.omega,g.m=-wf,-mf
        if spin==1:
            ed=g.derivative(kernel,-2,0,'eth')[0]+g.sc['tau']*kernel
            th=g.derivative(kernel,-2,0,'thorn')[0]+g.sc['rho']*kernel
            vec=[(-g.cov[0][i]*ed+g.cov[2][i]*th)/2 for i in range(4)]
        else:
            ed=g.derivative(kernel,2,0,'ethp')[0]+g.sc['taup']*kernel
            th=g.derivative(kernel,2,0,'thornp')[0]+g.sc['rhop']*kernel
            vec=[(g.cov[1][i]*ed-g.cov[3][i]*th)/2 for i in range(4)]
        if chirality<0:
            vec=[v.conjugate() for v in vec]
        g.omega,g.m=-omega,-m
        se=current_adjoint_contraction(g,vec,u,fup)
        dkw=dkw_current_adjoint(g,vec)
        correction=4j/(9*omega)*sum(u[i]*u[j]*dkw[i][j] for i in range(4) for j in range(4))
        results.append((8*np.pi*se.value/(ut*w0),8*np.pi*correction.value/(ut*w0)))
    return tuple(results)
