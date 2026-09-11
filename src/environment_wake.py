"""Assemble the complex scalar perturbation with its actual spheroidal modes."""
import numpy as np
from environment_response import SampledResponse
from environment_source import angular_mode


class EnvironmentalWake:
    def __init__(self,reports,ellmin=2,ellmax=6,allow_partial=False,
                 allow_mixed_discretization=False):
        if not 2<=ellmin<=ellmax:raise ValueError('This wake excludes the resonant cloud ell=1')
        self.modes={};reference=None;numerical_reference=None
        self.mode_provenance=[];self.mixed_discretization=False
        for data in reports:
            p=data['parameters'];metric=p['metric'];key=(p['scalar_ell'],p['scalar_m'])
            if not ellmin<=key[0]<=ellmax:continue
            if p.get('background') is not None:raise ValueError('This assembler requires a stationary Kerr cloud')
            if key in self.modes:raise ValueError(f'Duplicate scalar mode {key}')
            common=dict(alpha=p['alpha'],a=metric['a'],r0=metric['orbital_radius'],
                        cloud_mass=p['cloud_mass'],coordinates=metric['coordinates'],
                        flux_scaling=data.get('flux_scaling'))
            numerical=dict(metric_ellmax=metric['ellmax'],
                        angular_backend=metric.get('angular_backend','pybhpt-default'),
                        angular_order=p['angular_order'],radial_order=p['radial_order'],
                        source_panels=p['source_panels'],green_outer_radius=p['green_outer_radius'],
                        green_horizon_offset=p['green_horizon_offset'],
                        horizon_order=p.get('horizon_quadrature_order'),
                        horizon_log=p.get('horizon_log_first_panel',False),
                        infinity_method=p.get('infinity_method','series'))
            if reference is not None and common!=reference:
                raise ValueError('Mode reports have different physics or normalization')
            if numerical_reference is not None and numerical!=numerical_reference:
                if not allow_mixed_discretization:
                    raise ValueError('Mode reports have different discretization; explicit opt-in required')
                self.mixed_discretization=True
            reference=common
            if numerical_reference is None:numerical_reference=numerical
            if metric['m_g']!=key[1]-1 or (sum(key)%2):
                raise ValueError('Mode is inconsistent with the reflection-even |211> cloud')
            self.modes[key]=(p['omega'],SampledResponse.from_report(data))
            self.mode_provenance.append(dict(ell=key[0],m=key[1],numerical=numerical,
                                             metric=metric))
        if not self.modes:raise ValueError('No modes in the requested range')
        self.parameters=dict(reference)
        for key,value in numerical_reference.items():
            self.parameters[key]=(value if all(row['numerical'][key]==value
                for row in self.mode_provenance) else None)
        self.radial_domain=(max(response.green.rmin for _,response in self.modes.values()),
                            min(response.green.rmax for _,response in self.modes.values()))
        if self.radial_domain[0]>=self.radial_domain[1]:
            raise ValueError('No common radial solution domain')
        self.omega_p=1/(reference['r0']**1.5+reference['a'])
        omegas=[omega-(m-1)*self.omega_p for (ell,m),(omega,response) in self.modes.items()]
        if np.ptp(omegas)>1e-13:raise ValueError('Inconsistent background cloud frequency')
        self.omega_c=omegas[0]
        expected={(ell,m) for ell in range(ellmin,ellmax+1) for m in range(-ell,ell+1) if (ell+m)%2==0}
        self.missing=sorted(expected-set(self.modes))
        if self.missing and not allow_partial:
            raise ValueError(f'Missing scalar modes for full wake: {self.missing}')
        self.status='partial_mode_field' if self.missing else 'complete_mode_range_not_converged'

    def evaluate(self,r,theta,phi,time=0.):
        r,theta,phi=np.broadcast_arrays(np.asarray(r,float),np.asarray(theta,float),np.asarray(phi,float))
        if np.any(theta<0) or np.any(theta>np.pi):raise ValueError('Require polar angle in [0,pi]')
        if hasattr(self,'radial_domain') and (np.any(r<self.radial_domain[0]) or np.any(r>self.radial_domain[1])):
            raise ValueError('Requested radii outside the common radial solution domain')
        unique,index=np.unique(r.ravel(),return_inverse=True)
        result=np.zeros(r.shape,complex)
        for (ell,m),(omega,response) in self.modes.items():
            radial=response.evaluate(unique)[0][index].reshape(r.shape)
            # angular_mode also calculates a derivative with coordinate-pole
            # 0/0; its returned value is regular and has the exact axis limit.
            with np.errstate(divide='ignore',invalid='ignore'):
                angular=angular_mode(theta,ell,m,self.parameters['a']**2*(omega**2-self.parameters['alpha']**2))[0]
            if not np.all(np.isfinite(angular)):raise ValueError('Non-finite spheroidal value')
            result+=radial*angular*np.exp(1j*(m*phi-omega*time))
        return result
