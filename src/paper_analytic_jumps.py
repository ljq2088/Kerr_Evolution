"""2023 analytic spin-2 particle jumps, with no radial ODE solve.

Use the positive-spin source formulas with m -> -m for P_minus2.
The separation called Lambda in the paper equals lambda_plus2 + 4.
"""
from functools import lru_cache
import mpmath as mp
import numpy as np
from paper_precise_angular import PreciseRealHarmonic


@lru_cache(maxsize=128)
def curvature_jumps(r0,a,m,ell):
    # Match the same physical double input, while retaining the arithmetic
    # through the particle-source calculation at extended precision.
    r=np.longdouble(r0);a=np.longdouble(a);op=1/(r**np.longdouble('1.5')+a)
    mm=-m;w=mm*op;de=r*r-2*r+a*a;dep=2*(r-1)
    gtt=-1+2/r;gtp=-2*a/r;gpp=r*r+a*a+2*a*a/r
    ut=1/np.sqrt(-gtt-2*op*gtp-op*op*gpp)
    E=-(gtt+op*gtp)*ut;L=(gtp+op*gpp)*ut
    aa=(E*(r*r+a*a)-a*L)/de;bb=L-a*E
    harmonic=PreciseRealHarmonic(2,ell,mm,a*w)
    with mp.workdps(50):theta=np.longdouble(str(mp.pi/2));pi=np.longdouble(str(mp.pi))
    S=harmonic(theta);dS=harmonic(theta,deriv=1);lam=harmonic.eigenvalue+4
    Q=mm*(1-a*op);K=(r*r+a*a)*w-a*mm;W=K/de;Wp=(2*r*w*de-K*dep)/de**2
    C2=bb*bb*S
    C1=2*bb*(1j*aa*(dS+Q*S)+bb*(-1j*W+1/r)*S)
    C0=2*aa*(1j*bb*(2/r-1j*W)+aa*(-Q+1j*a/r))*(dS+Q*S)+(aa*aa*lam+bb*bb*(1j*Wp-2j*W/r-W*W))*S
    Vminus4=de*(W*W+1j*Wp)-1j*dep*W+6j*w*r-lam
    # The paper's m reversal holds its angular phase fixed. In the
    # Condon--Shortley basis used here, the spin flip/m reversal supplies
    # (-1)^(m+1), including the sign of the rescaled negative Weyl scalar.
    pref=(-1)**(m+1)*(-4*pi/(r*r*ut)*de)
    return lam,pref*(C1-dep/de*C2),pref*(C0-Vminus4/de*C2)
