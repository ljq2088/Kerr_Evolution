"""Development bridge from independently reconstructed Lorenz modes to the cloud.

Only nonstatic modes are available here. Static completion,
full source matching, and production convergence remain external requirements.
"""
from functools import lru_cache
import numpy as np
from lorenz_metric import nonstatic_metric


class LorenzMetricMode:
    def __init__(self,orbital_radius,a,m,ellmax):
        if abs(m)<1 or ellmax<abs(m):
            raise ValueError('This development adapter requires nonzero m_g and ellmax>=|m_g|')
        self.r0,self.a,self.m,self.ellmax=orbital_radius,a,m,ellmax
        self.omega=m/(orbital_radius**1.5+a)

    @lru_cache(maxsize=8192)
    def _values(self,r,theta):
        h=np.zeros((4,4),complex)
        for ell in range(abs(self.m),self.ellmax+1):
            _,piece=nonstatic_metric(r,theta,self.r0,self.a,ell,self.m,order=6)
            h+=np.array([[v.value for v in row] for row in piece])
        return tuple(h.ravel())

    def __call__(self,r,theta):
        if r==self.r0:
            raise ValueError('Use one-sided source sampling at the orbit')
        return np.asarray(self._values(float(r),float(theta))).reshape(4,4).copy()

    @property
    def provenance(self):
        result=dict(status='development_not_production_validated',coordinates='Boyer-Lindquist',
                    field='covariant h_ab; particle mass stripped',m_g=self.m,ellmax=self.ellmax,
                    orbital_radius=self.r0,a=self.a,static_and_low_modes_included=False,
                    kappa_mass_squared_step=min(5e-5,.01*self.omega**2),
                    maxwell_chiralities='both; full complex compact current')
        if abs(self.m)==1:
            result['nonstatic_dipole_included']=True
        return result
