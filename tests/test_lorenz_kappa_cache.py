"""Prevent cyclic L=20 trace/kappa requests from evicting all radial solutions."""
from types import SimpleNamespace
import numpy as np
import pytest
import lorenz_kappa as kappa

def test_full_L20_angular_sweep_reuses_radial_solutions(monkeypatch):
    builds=[]
    def fake_green(a,mu,omega,ell,m,**kwargs):
        builds.append((ell,kwargs['mass_squared']))
        value=np.array([complex(ell,kwargs['mass_squared']),1.+0j])
        return SimpleNamespace(insol=SimpleNamespace(sol=lambda r:value.copy()),
            upsol=SimpleNamespace(sol=lambda r:2*value.copy()),w0=3.+0j)
    monkeypatch.setattr(kappa,'RadialGreen',fake_green)
    monkeypatch.setattr(kappa,'angular_mode',lambda *args:(1.,0.,2.))
    kappa.scalar_resolvent.cache_clear()
    try:
        omega=1/(20**1.5+.88);step=min(5e-5,.01*omega**2)
        first={}
        for radius in (3.,10.):
            for theta_index in range(40):
                for ell in range(1,21):
                    for msq in (0.,step,-step,step/2,-step/2):
                        value=kappa.scalar_resolvent(20.,.88,ell,1,msq,1e-12,2000.,'series')
                        key=(ell,msq)
                        if key in first:assert value is first[key]
                        else:first[key]=value
        info=kappa.scalar_resolvent.cache_info()
        assert len(builds)==100
        assert info.misses==100 and info.hits==7900
    finally:kappa.scalar_resolvent.cache_clear()

@pytest.mark.parametrize('msq',[0.,1e-6,-1e-6])
def test_cached_solution_equals_fresh_radial_solve(msq):
    kappa.scalar_resolvent.cache_clear()
    try:
        args=(20.,.88,2,1,msq,1e-12,2000.,'series')
        cached=kappa.scalar_resolvent(*args)
        fresh=kappa.scalar_resolvent.__wrapped__(*args)
        assert cached[1:]==fresh[1:]
        radii=np.array([2.,10.,20.,100.,320.])
        for name in ('insol','upsol'):
            np.testing.assert_array_equal(getattr(cached[0],name).sol(radii),
                                          getattr(fresh[0],name).sol(radii))
    finally:kappa.scalar_resolvent.cache_clear()
