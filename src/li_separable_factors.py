"""Li 2507.02045v2 separable Gamma factors (Fourier p<=12, 64 digits).

This is an independent building block, not the full environmental source.
Angular coefficients can be precomputed independently of radial factors.
No old numerical-grid source is relabelled as a semi-analytical Li source.
"""
from functools import lru_cache
from math import comb
import mpmath as mp


def _integer(value, name, lower=0, upper=None):
    if isinstance(value,bool) or int(value)!=value or value<lower or (upper is not None and value>upper):
        raise ValueError(f'Invalid {name}')
    return int(value)


@lru_cache(maxsize=8192)
def _coefficients(r_text,a_text,beta,pmax,dps):
    with mp.workdps(dps):
        r,a=mp.mpf(r_text),mp.mpf(a_text)
        if not mp.isfinite(r) or not mp.isfinite(a) or r<=0:
            raise ValueError('Finite r>0 and finite a required')
        A=r*r+a*a/2;B=a*a/2
        if beta==0:return tuple(mp.mpf(int(p==0)) for p in range(pmax+1))
        def base(A,n):
            D=mp.sqrt(A*A-B*B)
            t=-B/(A+D)
            return (1 if n==0 else 2)*t**n/D
        return tuple(mp.mpf(0) if p%2 else
            (-1)**(beta-1)*mp.diff(lambda z:base(z,p//2),A,beta-1)/mp.factorial(beta-1)
            for p in range(pmax+1))


def sigma_fourier_coefficients(r,a,beta,pmax=12,dps=64):
    """Return f_p for (r^2+a^2 cos(theta)^2)^(-beta).

    Precision applies to arithmetic on supplied inputs; decimal strings retain
    more input digits than binary floats. It does not upgrade input accuracy.
    Beta0..3 is the range needed by the published source construction.
    """
    beta=_integer(beta,'beta',0,3);pmax=_integer(pmax,'pmax');dps=_integer(dps,'dps',20)
    return _coefficients(str(r),str(a),beta,pmax,dps)


def gamma_product_terms(r,a,beta,sigma,pmax=12,dps=64):
    """Terms (p,k,c) with Gamma^-beta barGamma^-sigma = sum c cos(p*t) cos(t)^k.

    Equality is up to Sigma Fourier truncation. Includes binomial coefficients,
    which are missing from the displayed binomial identity in the paper.
    Gamma=r+i a cos(theta). All complex coefficients are preserved.
    """
    beta=_integer(beta,'beta',0,3);sigma=_integer(sigma,'sigma',0,3)
    dps=_integer(dps,'dps',20)
    f=sigma_fourier_coefficients(r,a,max(beta,sigma),pmax,dps)
    with mp.workdps(dps):
        r,a=mp.mpf(str(r)),mp.mpf(str(a));n=abs(beta-sigma)
        sign=-1 if beta>=sigma else 1
        return tuple((p,k,fp*comb(n,k)*r**(n-k)*(sign*1j*a)**k)
                     for p,fp in enumerate(f) if fp for k in range(n+1))


def evaluate_terms(terms,theta,dps=64):
    with mp.workdps(dps):
        t=mp.mpf(str(theta));x=mp.cos(t)
        return mp.fsum(c*mp.cos(p*t)*x**k for p,k,c in terms)
