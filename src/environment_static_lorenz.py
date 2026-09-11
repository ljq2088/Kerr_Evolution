"""Static Lorenz field candidate assembled from independently matched jumps.

Auxiliary fields use regular In/Up Green solutions. Completion conventions
are explicit: Berndtson has shifted interior charges; physical fixes them to
zero but is not asymptotically flat in these Lorenz coordinates. Finite-mode
source matching and boundary convergence must still be checked per dataset.
"""
import json
from functools import lru_cache
from pathlib import Path
import numpy as np
from lorenz_static_spin2 import sourced_static_spin2
from lorenz_static_gauge import static_trace_metric
from lorenz_static_boundary import matched_free_scalar_metric,regular_completion_y2
from lorenz_completion import completion_metric


class StaticLorenzMode:
    def __init__(self,matching,completion='berndtson'):
        if isinstance(matching,(str,Path)):
            matching=json.loads(Path(matching).read_text())
        if completion not in ('berndtson','physical'):
            raise ValueError('Completion must be berndtson or physical')
        meta=matching['sample_metadata']
        if not meta['circular_isometry_average'] or matching['charge_conditions']!='exact_elimination':
            raise ValueError('Require circular matching with exact conserved charges')
        self.a,self.r0=meta['a'],meta['r0']
        self.ellmax=matching['parameters']['ellmax']
        self.m,self.omega,self.ellmin=0,0.,0
        self.completion=completion
        self.jumps=dict(zip(matching['basis_labels'],matching['coefficients']))
        # The matching used y2(r0)=y2'(r0)=0. Removing its growing homogeneous
        # mode changes D/E; compensate the free scalar jump in unit Y20.
        y2=regular_completion_y2(self.a,self.r0)
        norm=np.sqrt(5/(4*np.pi))
        for datum in (0,1):
            self.jumps[f'kappa_2_{datum}']-=2*(self.jumps['D']-self.jumps['E'])*float(y2[datum].real)/norm
        b=np.sqrt(1-self.a*self.a);rp,rm=1+b,1-b
        inside={key:0. for key in 'BCDEFG'}
        inside['B']=-self.jumps['B']
        if completion=='berndtson':
            inside['D']=-self.jumps['D'];inside['F']=-self.jumps['F']
            inside['E']=-self.jumps['D']
            inside['G']=-self.a*self.jumps['D']-rp*(rp-rm)**2/(rp+rm)*self.jumps['F']
            # Fix horizon position as well as angular velocity in this
            # explicit completion basis. The printed polynomial coefficient
            # of Delta c_F does not cancel the horizon pole with the stated z.
            # C_F=-2*a*b*rp**2 follows from delta r_H=0; see STATIC_BOUNDARY.md.
            inside['C']=-2*rp*rp*self.jumps['D']+2*self.a*b*rp*rp*self.jumps['F']
        self.inside=inside
        self.outside={key:inside[key]+self.jumps[key] for key in inside}
        self.scalar_degrees=sorted(int(key.split('_')[1]) for key in self.jumps if key.startswith('kappa_') and key.endswith('_0'))

    def metric_jet(self,r,theta,order=8):
        if r==self.r0:
            raise ValueError('Select a one-sided vacuum radius')
        g,total=sourced_static_spin2(r,theta,r0=self.r0,a=self.a,ell=2,order=order,boundary='regular')
        def add(h,coefficient=1.):
            for i in range(4):
                for j in range(4):total[i][j]+=coefficient*h[i][j]
        for ell in range(2,self.ellmax+1):
            if ell>2:
                _,h=sourced_static_spin2(r,theta,r0=self.r0,a=self.a,ell=ell,order=order,boundary='regular')
                add(h)
            if ell%2==0:
                _,h,_=static_trace_metric(r,theta,r0=self.r0,a=self.a,ell=ell,order=order,boundary='regular')
                add(h)
        for ell in self.scalar_degrees:
            _,h,_=matched_free_scalar_metric(r,theta,self.r0,self.a,ell,
                self.jumps[f'kappa_{ell}_0'],self.jumps[f'kappa_{ell}_1'],order=order)
            add(h)
        for key,value in (self.inside if r<self.r0 else self.outside).items():
            if value:
                _,h,_=completion_metric(r,theta,a=self.a,mode=key,order=order,
                                       reference_radius=self.r0,y_boundary='regular-quadrupole')
                add(h,value)
        return g,total

    @lru_cache(maxsize=8192)
    def _values(self,r,theta):
        _,h=self.metric_jet(r,theta,order=6)
        return tuple(v.value for row in h for v in row)

    def __call__(self,r,theta):
        return np.asarray(self._values(float(r),float(theta))).reshape(4,4).copy()

    @property
    def provenance(self):
        return dict(status='development_static_boundary_candidate',coordinates='Boyer-Lindquist',
                    m_g=0,a=self.a,orbital_radius=self.r0,ellmax=self.ellmax,
                    completion=self.completion,auxiliary_boundary='regular In/Up',
                    y_quadrupole='regular In/Up; free scalar jumps transformed',
                    matched_jumps=self.jumps,
                    inside_completion=self.inside,outside_completion=self.outside,
                    full_source_convergence_verified=False)
