"""Weyl inputs for Lorenz reconstruction, sourced via pinned pybhpt.

These are curvature fields, not a claim of a completed Lorenz metric.
"""
from pybhpt.geo import KerrGeodesic
from pybhpt.teuk import TeukolskyMode
from functools import lru_cache


@lru_cache(maxsize=256)
def weyl_amplitudes(r0, a=.6, ell=2, m=2):
    geo=KerrGeodesic(a,r0,0.,1.,nsamples=128)
    result={}
    for s in (-2,2):
        mode=TeukolskyMode(s,ell,m,0,0,geo)
        mode.solve(geo)
        result[s]=(mode.amplitude('Up'),mode.amplitude('In'))
    return result


def integrated_weyl_amplitudes(r0, a=.6, ell=2, m=2):
    """Nonstatic time integrals, retained separately from physical curvature.

    With exp(-i omega t), inverse Lie_T is i/omega. The public 2406.12510v3
    amplitude table agrees with these time integrals rather than raw pybhpt
    curvature. This empirical convention bridge must remain explicit.
    """
    omega=m/(r0**1.5+a)
    if omega==0:
        raise ValueError('Static modes require separate completion; no inverse time derivative')
    return {s:tuple(1j*z/omega for z in pair)
            for s,pair in weyl_amplitudes(r0,a,ell,m).items()}


if __name__=='__main__':
    for r in (4.,6.):
        print(r,weyl_amplitudes(r))
