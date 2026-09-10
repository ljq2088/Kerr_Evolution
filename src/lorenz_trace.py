"""First of the six Lorenz reconstruction fields: trace from a circular particle.

Box h = 16 pi T, T = -mp delta(r-r0) delta(theta-pi/2)
delta(phi-Omega t)/(ut Sigma sin theta). M=mp=1 (q stripped).
This trace alone is not the Lorenz metric and must not be used as one.
"""
import numpy as np
from environment_source import angular_mode, kerr_metric
from environment_radial import RadialGreen


def trace_amplitudes(r0, a=.6, ell=2, m=2, rmax=4000., offset=1e-6):
    op = 1/(r0**1.5+a)
    omega = m*op
    g = kerr_metric(r0,np.pi/2,a)
    ut = 1/np.sqrt(-g[0,0]-2*op*g[0,3]-op*op*g[3,3])
    s = angular_mode(np.pi/2,ell,m,a*a*omega*omega)[0]
    J = -16*np.pi*s/ut
    radial = RadialGreen(a,0.,omega,ell,m,rmax=rmax,offset=offset,rtol=1e-11)
    zh = radial.upsol.sol(r0)[0]*J/radial.w0
    zi = radial.insol.sol(r0)[0]*J/radial.w0
    # Keep the same tortoise phase convention as the published amplitude table.
    return zi, zh


if __name__ == '__main__':
    for r in (4.,6.,10.,20.):
        zi,zh=trace_amplitudes(r)
        print(r,zi,zh,abs(zi),abs(zh))
