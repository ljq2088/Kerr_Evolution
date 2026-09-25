"""Li-style separated tetrad source projection with cached angular couplings.

Consumes weighted SPHERICAL metric coefficients, never bare spheroidal modes.
The Appendix's proven missing term is included by default; printed_only is an
explicit diagnostic. No claim that the authors used the erroneous printed term.
Temporal/cloud frequency treatment and mass normalization belong to the caller.
"""
from functools import lru_cache
from fractions import Fraction
from pathlib import Path
import json
import numpy as np
from pybhpt.swsh import Yslm
from li_separable_factors import gamma_product_terms

SPINS=(0,0,2,-2,1,-1,1,-1,0,0)
TERMS=json.loads(Path(__file__).with_name('li_source_terms.json').read_text())

class LiSeparatedSource:
    def __init__(self,*,a,omega_c,m_c,metric_m,spherical_lmax,cloud_angular,target_angular,
                 quadrature=64,pmax=12,dps=64,printed_only=False):
        if quadrature<2 or spherical_lmax<abs(metric_m):raise ValueError('Invalid angular grid/truncation')
        self.a,self.w,self.mc=a,omega_c,m_c
        self.mg,self.lmax,self.pmax,self.dps=metric_m,spherical_lmax,pmax,dps
        x,self.weights=np.polynomial.legendre.leggauss(quadrature)
        self.theta=np.arccos(x);self.s=np.sin(self.theta);self.c=x
        self.cloud=np.asarray(cloud_angular(self.theta),complex)
        self.target=np.asarray(target_angular(self.theta),complex)
        if self.cloud.shape!=(3,quadrature) or self.target.shape!=(quadrature,):
            raise ValueError('Expected three cloud angular derivatives and target conjugate/dual weights')
        if not np.isfinite(self.cloud).all() or not np.isfinite(self.target).all():raise ValueError('Nonfinite angular input')
        self.ys={(j,c):np.asarray(Yslm(spin,j,metric_m,self.theta),complex)
                 for c,spin in enumerate(SPINS) for j in range(max(abs(metric_m),abs(spin)),spherical_lmax+1)}
        self.terms=TERMS['printed_terms']+([] if printed_only else TERMS['required_correction_terms'])
        self.printed_only=printed_only

    @lru_cache(maxsize=None)
    def _angular(self,component,j,derivative,sin_power,cos_power,p):
        values=self.target*self.cloud[derivative]*self.ys[j,component]
        return 2*np.pi*np.dot(self.weights,values*self.s**sin_power*self.c**cos_power*np.cos(p*self.theta))

    def project(self,r,cloud_radial,metric_spherical):
        """Return J_lm(r)=integral target*Sigma*S dOmega at one radius.

        cloud_radial=[R,R',R'']; metric[j,component] follows the ten-component
        weighted tetrad schema, with trace already converted to spherical.
        """
        metric=np.asarray(metric_spherical,complex);R=np.asarray(cloud_radial,complex)
        if metric.shape!=(self.lmax+1,10) or R.shape!=(3,):raise ValueError('Wrong input shape')
        if not np.isfinite(metric).all() or not np.isfinite(R).all():raise ValueError('Nonfinite radial input')
        d=r*r-2*r+self.a**2;dp=2*(r-1)
        if r<=1+np.sqrt(1-self.a**2):raise ValueError('Exterior r required')
        K=(r*r+self.a**2)*self.w-self.a*self.mc;V=K/d
        Vp=(2*r*self.w*d-K*dp)/d**2
        radial_values=(r,self.a,d,dp,V,Vp,self.w,self.mc)
        factors={};result=0j
        for t in self.terms:
            key=t['beta'],t['sigma']
            if key not in factors:
                factors[key]=[(p,k,complex(c)) for p,k,c in gamma_product_terms(
                    r,self.a,*key,pmax=self.pmax,dps=self.dps)]
            coefficient=complex(float(Fraction(t['coefficient'][0])),float(Fraction(t['coefficient'][1])))
            coefficient*=R[t['radial_derivative']]
            for v,power in zip(radial_values,t['radial_powers']):
                if power:coefficient*=v**power
            comp=t['component']
            for j in range(max(abs(self.mg),abs(SPINS[comp])),self.lmax+1):
                if metric[j,comp]==0:continue
                coupling=sum(c*self._angular(comp,j,t['angular_derivative'],
                    t['sin_power'],t['cos_power']+k,p) for p,k,c in factors[key])
                result+=coefficient*metric[j,comp]*coupling
        return result
