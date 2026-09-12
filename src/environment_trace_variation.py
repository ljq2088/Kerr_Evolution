"""Analytic auxiliary-mass variation of the forced scalar trace resolvent.

Diagnostic implementation; not imported by production metric reconstruction.
Includes the varying angular source projector as well as both radial boundaries.
"""
import numpy as np
from environment_radial_variation import RadialMassVariation
from environment_angular_variation import angular_mode_mass_derivative
from environment_source import kerr_metric


class TraceMassVariation:
    def __init__(self,r0,a,ell,m,*,mass_squared=0.,rmax=2000.,rtol=1e-12):
        self.r0=r0;self.a=a;self.ell=ell;self.m=m;self.mass_squared=mass_squared
        self.omega=m/(r0**1.5+a)
        self.radial=RadialMassVariation(a,self.omega,ell,m,mass_squared=mass_squared,
                                      rmax=rmax,rtol=rtol)
        if not self.radial.rmin<r0<rmax:
            raise ValueError('Particle must be inside radial integration domain')
        op=1/(r0**1.5+a)
        metric=kerr_metric(r0,np.pi/2,a)
        ut=1/np.sqrt(-metric[0,0]-2*op*metric[0,3]-op*op*metric[3,3])
        s,_,_,ds,_,_=angular_mode_mass_derivative(np.pi/2,ell,m,a,self.omega,mass_squared)
        self.source=-16*np.pi*s/ut;self.source_mass=-16*np.pi*ds/ut
        w,dw=self.radial.w0,self.radial.w0_mass
        def amplitude(solution):
            r,_,dr,_=solution.sol(r0)
            z=self.source*r/w
            dz=(self.source_mass*r+self.source*dr)/w-z*dw/w
            return z,dz
        self.zi,self.zi_mass=amplitude(self.radial.insol)
        self.zh,self.zh_mass=amplitude(self.radial.upsol)

    def radial_state(self,r):
        """Return F,F_r,F_nu,F_rnu; choose the exterior derivative at r=r0."""
        r=np.asarray(r,float)
        if np.any(r<self.radial.rmin) or np.any(r>self.radial.rmax):
            raise ValueError('Evaluation outside integrated radial domain')
        inner=self.radial.insol.sol(r);outer=self.radial.upsol.sol(r)
        z=np.where(r<self.r0,self.zh,self.zi)
        dz=np.where(r<self.r0,self.zh_mass,self.zi_mass)
        state=np.where(r<self.r0,inner,outer)
        return np.array([z*state[0],z*state[1],
                         dz*state[0]+z*state[2],dz*state[1]+z*state[3]])

    def field_state(self,r,theta):
        """h,h_r,h_theta and their nu derivatives, without exp(-iwt+imphi)."""
        f,fr,df,dfr=self.radial_state(r)
        s,st,_,ds,dst,_=angular_mode_mass_derivative(
            theta,self.ell,self.m,self.a,self.omega,self.mass_squared)
        return np.array([f*s,fr*s,f*st,df*s+f*ds,dfr*s+fr*ds,df*st+f*dst])

    def field_jet(self,g):
        """Return trace jet and its mass derivative by differentiated ODE recurrences."""
        from lorenz_jet import Jet
        if g.a!=self.a or g.m!=self.m or not np.isclose(g.omega,self.omega,rtol=1e-14,atol=0):
            raise ValueError('Geometry frequency and mode must match trace resolvent')
        n=g.order;r=g.r;t=g.theta
        f,fr,df,dfr=self.radial_state(float(r.value.real))
        s,st,A,ds,dst,dA=angular_mode_mass_derivative(
            float(t.value.real),self.ell,self.m,self.a,self.omega,self.mass_squared)
        R=Jet(f,n);D=Jet(df,n);S=Jet(s,n);E=Jet(ds,n)
        R.c[1,0]=fr;D.c[1,0]=dfr;S.c[0,1]=st;E.c[0,1]=dst
        K=(r*r+self.a*self.a)*self.omega-self.a*self.m
        V=K*K/g.delta-self.mass_squared*r*r-self.a**2*self.omega**2+2*self.a*self.m*self.omega-A
        U=self.a**2*(self.omega**2-self.mass_squared)*t.cos()**2+A-self.m**2/t.sin()**2
        for j in range(n-1):
            factor=(j+1)*(j+2)
            rhs=-(2*(r-1)*R.derivative(0)+V*R)/g.delta
            drhs=(-2*(r-1)*D.derivative(0)-V*D+(r*r+dA)*R)/g.delta
            R.c[j+2,0]=rhs.c[j,0]/factor
            D.c[j+2,0]=drhs.c[j,0]/factor
            rhs=-t.cos()/t.sin()*S.derivative(1)-U*S
            drhs=-t.cos()/t.sin()*E.derivative(1)-U*E+(self.a**2*t.cos()**2-dA)*S
            S.c[0,j+2]=rhs.c[0,j]/factor
            E.c[0,j+2]=drhs.c[0,j]/factor
        return R*S,D*S+R*E

    def kappa_jet(self,g):
        return -1j*self.omega*self.field_jet(g)[1]


def install_analytic_kappa_diagnostic():
    """Process-local replacement for comparison runs; leaves source files unchanged."""
    from functools import lru_cache
    import lorenz_metric
    @lru_cache(maxsize=96)
    def solution(r0,a,ell,m,rtol,rmax):
        return TraceMassVariation(r0,a,ell,m,rtol=rtol,rmax=rmax)
    def analytic(g,r0,ell,step=None,rtol=1e-12,rmax=2000.,infinity_method='series'):
        if step is not None or infinity_method!='series':
            raise ValueError('Analytic diagnostic uses series boundaries and no difference step')
        return solution(r0,g.a,ell,g.m,rtol,rmax).kappa_jet(g)
    lorenz_metric.kappa_jet=analytic
