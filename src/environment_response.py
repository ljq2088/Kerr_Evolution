"""Continuous response from completed, piecewise Gauss source samples.

Interpolate J on each original panel, not the resulting field. Integrate the
two Green amplitudes from opposite ends to avoid subtracting nearly equal
outer integrals. No source is invented outside the recorded finite cutoffs.
"""
import numpy as np
from scipy.integrate import solve_ivp
from numpy.polynomial.legendre import legfit, legval
from environment_radial import RadialGreen


class SampledResponse:
    def __init__(self, green, panels, radii, source, *, log_first=False,
                 rtol=1e-10, atol=1e-14, phase=None):
        self.green=green
        self.panels=np.asarray(panels,float)
        r=np.asarray(radii,float);J=np.asarray(source,complex)
        if (len(r)!=len(J) or np.any(np.diff(r)<=0)
            or np.any(np.diff(self.panels)<=0)
            or self.panels[0]<green.rmin or self.panels[-1]>green.rmax
            or not np.all(np.isfinite(J))):
            raise ValueError('Invalid completed source grid')
        self.parts=[]
        count=0
        for index,(lo,hi) in enumerate(zip(self.panels[:-1],self.panels[1:])):
            mask=(r>lo)&(r<hi);nodes=r[mask];values=J[mask]
            count+=len(nodes)
            if len(nodes)<2:raise ValueError('Every panel needs at least two interior samples')
            logarithmic=bool(index==0 and log_first)
            z0,z1=(np.log([lo-green.rp,hi-green.rp]) if logarithmic else (lo,hi))
            z=np.log(nodes-green.rp) if logarithmic else nodes
            x=2*(z-z0)/(z1-z0)-1
            expected=np.polynomial.legendre.leggauss(len(nodes))[0]
            if not np.allclose(x,expected,rtol=0,atol=2e-9):
                raise ValueError('Samples do not match the declared Gauss panel')
            reduced=values if phase is None else values*np.exp(-1j*phase(nodes))
            coefficients=legfit(x,reduced,len(nodes)-1)
            def rhs(z,y,branch=0,coeff=coefficients,a=z0,b=z1,log=logarithmic):
                radius=green.rp+np.exp(z) if log else z
                jacobian=np.exp(z) if log else 1.
                source_value=legval(2*(z-a)/(b-a)-1,coeff)
                if phase is not None:source_value*=np.exp(1j*phase(radius))
                solution=green.insol if branch==0 else green.upsol
                return [jacobian*solution.sol(radius)[0]*source_value/green.w0]
            # Bound channels can have tiny Green coefficients multiplying huge
            # homogeneous solutions. Scale each integral before imposing atol.
            jac=nodes-green.rp if logarithmic else np.ones(len(nodes))
            left_scale=float(np.max(abs(green.insol.sol(nodes)[0]*values*jac/green.w0))*(z1-z0)) or 1.
            right_scale=float(np.max(abs(green.upsol.sol(nodes)[0]*values*jac/green.w0))*(z1-z0)) or 1.
            left=solve_ivp(lambda z,y: rhs(z,y,0)[0]/left_scale,(z0,z1),[0j],method='DOP853',rtol=rtol,atol=atol,dense_output=True)
            right=solve_ivp(lambda z,y: -rhs(z,y,1)[0]/right_scale,(z1,z0),[0j],
                method='DOP853',rtol=rtol,atol=atol,dense_output=True)
            if not left.success or not right.success:
                raise RuntimeError('Cumulative source quadrature failed')
            self.parts.append(dict(left=left,right=right,log=logarithmic,left_scale=left_scale,right_scale=right_scale))
        if count!=len(r):raise ValueError('Source nodes outside open panel intervals')
        left_total=np.array([p['left'].y[0,-1]*p['left_scale'] for p in self.parts])
        right_total=np.array([p['right'].y[0,-1]*p['right_scale'] for p in self.parts])
        self.prefix=np.r_[0j,np.cumsum(left_total)]
        self.suffix=np.r_[np.cumsum(right_total[::-1])[::-1],0j]
        self.up_coefficient=self.prefix[-1]
        self.horizon_coefficient=self.suffix[0]

    def evaluate(self,radii):
        r=np.atleast_1d(radii).astype(float)
        if r.ndim!=1 or np.any(r<self.green.rmin) or np.any(r>self.green.rmax):
            raise ValueError('Requested radii outside the radial solution')
        A=np.zeros(r.size,complex);B=np.zeros(r.size,complex)
        A[r>=self.panels[-1]]=self.prefix[-1]
        B[r<self.panels[0]]=self.suffix[0]
        for index,p in enumerate(self.parts):
            mask=(r>=self.panels[index])&(r<self.panels[index+1])
            z=np.log(r[mask]-self.green.rp) if p['log'] else r[mask]
            if len(z):
                A[mask]=self.prefix[index]+p['left'].sol(z)[0]*p['left_scale']
                B[mask]=self.suffix[index+1]+p['right'].sol(z)[0]*p['right_scale']
        u,du=self.green.insol.sol(r);v,dv=self.green.upsol.sol(r)
        return np.asarray([v*A+u*B,dv*A+du*B])

    @classmethod
    def from_report(cls,data,dephase_outgoing=False,**kwargs):
        if data.get('status')!='truncated_single_mode_not_converged' or 'flux' not in data:
            raise ValueError('A completed finite-source response report is required')
        if data.get('flux_validity')=='historical_unreliable_boundary_result':
            raise ValueError('Historical invalid radial boundary cannot be used')
        p=data['parameters'];m=p['metric']
        green=RadialGreen(m['a'],p['alpha'],p['omega'],p['scalar_ell'],p['scalar_m'],
            rmax=p['green_outer_radius'],offset=p['green_horizon_offset'],rtol=1e-11,
            infinity_method=p.get('infinity_method','series'))
        if green.infinity_method=='series' and green.series_last_term_relative>1e-3:
            raise ValueError('Unresolved infinity series')
        rows=data['samples']
        phase=None
        if dephase_outgoing:
            from environment_radial import tortoise
            orbit=m['orbital_radius'];a=m['a'];frequency=m['m_g']/(orbit**1.5+a)
            origin=tortoise(orbit,a)
            # The particle is already a panel boundary. Match the carrier
            # continuously to zero there; inner source interpolation is kept.
            phase=lambda radius:frequency*(tortoise(np.maximum(radius,orbit),a)-origin)
        return cls(green,p['source_panels'],[s['r'] for s in rows],
            [complex(*s['source']) for s in rows],
            log_first=p.get('horizon_log_first_panel',False),phase=phase,**kwargs)
