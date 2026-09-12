"""Analytic scalar angular derivatives with respect to auxiliary mass squared.

Diagnostic support for a differentiated retarded scalar resolvent. This is
not yet used by production Lorenz reconstruction. Real parameters only.
"""
from functools import lru_cache
import numpy as np
from scipy.special import gammaln,lpmv
from environment_cloud import angular_basis_size


@lru_cache(maxsize=128)
def angular_eigenpair_mass_derivative(ell,m,a,omega,mass_squared=0.,size=20):
    """Return A, coefficients, dA/dmass2, dcoeff/dmass2 and harmonic degrees."""
    size=angular_basis_size(ell,m,size)
    if not all(np.isreal(v) and np.isfinite(v) for v in (a,omega,mass_squared)):
        raise ValueError('Real finite parameters required')
    ell,m=int(ell),abs(int(m))
    degrees=np.arange(m,m+size+1)
    cosine=np.zeros((size+1,size+1))
    for j,l in enumerate(degrees[:-1]):
        cosine[j,j+1]=cosine[j+1,j]=np.sqrt(((l+1)**2-m*m)/((2*l+1)*(2*l+3)))
    cosine2=(cosine@cosine)[:size,:size]
    derivative=a*a*cosine2
    matrix=np.diag(degrees[:size]*(degrees[:size]+1.))-a*a*(omega*omega-mass_squared)*cosine2
    eigenvalues,vectors=np.linalg.eigh(matrix)
    index=ell-m;eigenvalue=eigenvalues[index]
    coefficient=vectors[:,index].copy()
    coefficient*=1 if coefficient[index]>=0 else -1
    slope=float(coefficient@derivative@coefficient)
    gaps=eigenvalue-eigenvalues
    other=np.arange(size)!=index
    if np.any(abs(gaps[other])<1e-12):raise ValueError('Degenerate eigenpair derivative requires a subspace treatment')
    projected=vectors.T@derivative@coefficient
    amplitudes=np.zeros(size);amplitudes[other]=projected[other]/gaps[other]
    tangent=vectors@amplitudes
    # The normalization/phase convention is v.T vdot=0.
    tangent-=coefficient*(coefficient@tangent)
    for array in (coefficient,tangent):array.setflags(write=False)
    ls=degrees[:size].copy();ls.setflags(write=False)
    return float(eigenvalue),coefficient,slope,tangent,ls


def angular_mode_mass_derivative(theta,ell,m,a,omega,mass_squared=0.,size=20):
    """S, dtheta S, A, dmass2 S, dtheta dmass2 S, dmass2 A.

    Uses unit sphere norm including 2*pi in phi. Polar derivatives at the
    coordinate axes require separate limiting formulas and are rejected.
    """
    value,coefficient,slope,tangent,degrees=angular_eigenpair_mass_derivative(ell,m,a,omega,mass_squared,size)
    theta=np.asarray(theta,float)
    if not np.all(np.isfinite(theta)) or np.any(theta<=0) or np.any(theta>=np.pi):
        raise ValueError('Polar derivative evaluation requires 0 < theta < pi')
    mm=abs(int(m));phase=(-1)**mm if m<0 else 1
    x=np.cos(theta);st=np.sin(theta)
    outputs=[np.zeros_like(x) for _ in range(4)]
    for l,b,db in zip(degrees,coefficient,tangent):
        norm=np.sqrt((2*l+1)/(4*np.pi)*np.exp(gammaln(l-mm+1)-gammaln(l+mm+1)))
        leg=lpmv(mm,l,x);prev=lpmv(mm,l-1,x) if l>mm else np.zeros_like(x)
        y=phase*norm*leg;dy=phase*norm*(l*x*leg-(l+mm)*prev)/st
        outputs[0]+=b*y;outputs[1]+=b*dy;outputs[2]+=db*y;outputs[3]+=db*dy
    return outputs[0],outputs[1],value,outputs[2],outputs[3],slope
