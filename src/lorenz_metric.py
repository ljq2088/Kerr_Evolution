"""Experimental nonstatic complex Lorenz reconstruction in vacuum regions.

Includes spin-2, spin-1, trace, compact chi and resolvent kappa pieces.
Source matching and static completion remain unverified.
Do not use this intermediate assembly as a production environmental metric.
"""
import numpy as np
from functools import lru_cache
from pybhpt.radial import RadialTeukolsky
from pybhpt.swsh import SpinWeightedSpheroidalHarmonic
from lorenz_ghp import KerrGHP
from lorenz_mode_jet import separated_jet
from lorenz_corrector import Weighted
from lorenz_weyl import weyl_amplitudes
from lorenz_chi import sdag_tensor,direct_chi
from lorenz_tensor import tensor_divergence,vector_covariant_derivative
from lorenz_spin1 import spin1_amplitudes,cky_tensor
from lorenz_kappa import trace_field_jet,kappa_jet
from lorenz_chi import chi_amplitudes


@lru_cache(maxsize=8192)
def _homogeneous_radial_data(spin,ell,m,a,omega,r,bc):
    """Immutable radial data; independent of theta and source amplitude."""
    radial=RadialTeukolsky(spin,ell,m,a,omega,np.array([r]))
    radial.solve(bc=bc)
    return radial.eigenvalue,radial.radialsolution(bc,0),radial.radialderivative(bc,0)


def homogeneous_field_jet(g,spin,ell,amplitude,bc):
    r=float(g.r.value.real)
    theta=float(g.theta.value.real)
    eigenvalue,value,derivative=_homogeneous_radial_data(spin,ell,g.m,g.a,g.omega,r,bc)
    angular=SpinWeightedSpheroidalHarmonic(spin,ell,g.m,g.a*g.omega)
    R0=amplitude*value
    R1=amplitude*derivative
    return separated_jet(g,spin,eigenvalue,R0,R1,
                          angular(theta),angular(theta,deriv=1))


def h_operator(g,f,which,p):
    dk,ek=('thornp','ethp') if which==0 else ('thorn','eth')
    v,pp,qq=g.derivative(f,p,0,dk)
    v=g.zeta.conjugate()*v
    return g.derivative(v,pp,qq,ek)[0]


def spin2_metric(r,theta,r0,a=.6,ell=2,m=2,order=10,return_parts=False):
    if m==0 or r==r0:
        raise ValueError('Nonstatic vacuum points only')
    omega=m/(r0**1.5+a)
    g=KerrGHP(r,theta,a,omega=omega,m=m,order=order)
    bc,index=('In',1) if r<r0 else ('Up',0)
    amplitudes=weyl_amplitudes(r0,a,ell,m)
    psi0=homogeneous_field_jet(g,2,ell,amplitudes[2][index],bc)
    psi4=homogeneous_field_jet(g,-2,ell,amplitudes[-2][index],bc)
    z0,z4=g.zeta**4*psi0,g.zeta**4*psi4
    sd4,sd0=sdag_tensor(g,z0,4),sdag_tensor(g,z4,0)
    aab=[[4*(sd4[i][j]-sd0[i][j])/3 for j in range(4)] for i in range(4)]
    H4,H0=h_operator(g,z4,4,-4),h_operator(g,z0,0,4)
    l,n,mm,mb=g.cov
    htwo=[[(4j/(9*omega)/g.zeta**2)*(
        (l[i]*mm[j]-l[j]*mm[i])*H4/2-(mb[i]*n[j]-mb[j]*n[i])*H0/2)
        for j in range(4)] for i in range(4)]
    # Store the two-form in the first-index-divergence convention. This
    # orientation is fixed by div(xi_DKW)=0 and the Lorenz vacuum identity,
    # both independently tested; reversing it fails those identities.
    divH=tensor_divergence(g,htwo)
    chi=-(direct_chi(g,z4,4,-4)-direct_chi(g,z0,0,4))/(18*omega**2)
    if return_parts:
        return g,aab,[g.zeta**2*x for x in divH],chi
    xi=[g.zeta**2*divH[i]-g.partial(chi,i) for i in range(4)]
    dxi=vector_covariant_derivative(g,xi)
    h=[[(aab[i][j]-dxi[i][j]-dxi[j][i])*1j/omega for j in range(4)] for i in range(4)]
    return g,h


def spin1_metric(r,theta,r0,a=.6,ell=2,m=2,order=8,return_vector=False,full_current=False):
    """Vacuum spin-1 piece; optionally reconstruct the full complex current.

    full_current=False exposes the self-dual circularity diagnostic. The
    complete metric must include the independently sourced opposite chirality.
    """
    if m==0 or r==r0:
        raise ValueError('Nonstatic vacuum points only')
    omega=m/(r0**1.5+a)
    g=KerrGHP(r,theta,a,omega=omega,m=m,order=order)
    bc,index=('In',1) if r<r0 else ('Up',0)
    phi0=homogeneous_field_jet(g,1,ell,spin1_amplitudes(r0,a,ell,m,1)[index],bc)
    phi2=homogeneous_field_jet(g,-1,ell,spin1_amplitudes(r0,a,ell,m,-1)[index],bc)
    l,n,mm,mb=g.cov
    f=[[(mb[i]*n[j]-mb[j]*n[i])*phi0+(l[i]*mm[j]-l[j]*mm[i])*phi2
        for j in range(4)] for i in range(4)]
    if full_current:
        from lorenz_spin1_chiral import chiral_amplitudes
        reverse=KerrGHP(r,theta,a,omega=-omega,m=-m,order=order)
        anti=[]
        for s in (1,-1):
            se,dkw=chiral_amplitudes(r0,a,ell,m,s,-1)[index]
            anti.append(homogeneous_field_jet(reverse,s,ell,(se-dkw).conjugate(),bc).conjugate())
        f=[[f[i][j]+(mm[i]*n[j]-mm[j]*n[i])*anti[0]
             +(l[i]*mb[j]-l[j]*mb[i])*anti[1] for j in range(4)] for i in range(4)]
    ky=cky_tensor(g)
    htwo=[[sum((ky[i][k]*g.inv[k][c]*f[c][j]-ky[j][k]*g.inv[k][c]*f[c][i])/2
                for k in range(4) for c in range(4))*1j/omega for j in range(4)] for i in range(4)]
    # Factor 2 recovers BOTH input self-dual Maxwell scalars from F=d(xi).
    # For the full source, the opposite chirality is explicitly included
    # above; it cannot be replaced by a conjugation of a complex current.
    xi=[-2*x for x in tensor_divergence(g,htwo)]
    if return_vector:
        return g,xi
    dxi=vector_covariant_derivative(g,xi)
    h=[[-(dxi[i][j]+dxi[j][i])*1j/omega for j in range(4)] for i in range(4)]
    return g,h


def spin0_metric(r,theta,r0,a=.6,ell=2,m=2,order=8,kappa_step=None,return_vector=False):
    """Trace, compact chi and retarded kappa; no static completion."""
    if m==0 or r==r0:
        raise ValueError('Nonstatic vacuum points only')
    omega=m/(r0**1.5+a)
    g=KerrGHP(r,theta,a,omega=omega,m=m,order=order)
    htrace=trace_field_jet(g,r0,ell)
    kappa=kappa_jet(g,r0,ell,step=kappa_step)
    bc,index=('In',1) if r<r0 else ('Up',0)
    chi=homogeneous_field_jet(g,0,ell,chi_amplitudes(r0,a,ell,m)[index],bc)
    ky=cky_tensor(g)
    xi=[sum(ky[i][j]*g.inv[j][k]*g.partial(htrace,k)/2 for j in range(4) for k in range(4))
        +g.partial(kappa-chi,i) for i in range(4)]
    if return_vector:
        return g,xi
    dxi=vector_covariant_derivative(g,xi)
    h=[[-(dxi[i][j]+dxi[j][i])*1j/omega for j in range(4)] for i in range(4)]
    return g,h


def nonstatic_metric(r,theta,r0,a=.6,ell=2,m=2,order=10):
    """Complex development assembly; source matching is not yet validated.

    For all components, reliable Taylor derivatives extend only through
    order-6: the DKW scalar and its two gradients consume six orders.
    Use order>=8 for first/second derivatives of every component. The tt
    component consumes fewer radial orders because its gradients are temporal.

    AAB fields need not obey h_m=conj(h_-m): the physical real-part Fourier
    coefficient is (h_m+conj(h_-m))/2. That projection alone proves neither
    source normalization nor continuity; see arXiv:2406.12510v3 Sec. III A.
    """
    if m==0 or ell<max(1,abs(m)) or r==r0:
        raise ValueError('Nonstatic vacuum modes with ell>=|m| are required')
    if ell>=2:
        g,h2=spin2_metric(r,theta,r0,a,ell,m,order)
    else:
        # No spin-2 harmonic exists at ell=1. Its scalar and vector gauge
        # sectors must still be included in the nonstatic dipole mode.
        from lorenz_jet import Jet
        g=KerrGHP(r,theta,a,omega=m/(r0**1.5+a),m=m,order=order)
        h2=[[Jet(0.,order) for _ in range(4)] for _ in range(4)]
    _,h1=spin1_metric(r,theta,r0,a,ell,m,order,full_current=True)
    _,h0=spin0_metric(r,theta,r0,a,ell,m,order)
    return g,[[h2[i][j]+h1[i][j]+h0[i][j] for j in range(4)] for i in range(4)]
