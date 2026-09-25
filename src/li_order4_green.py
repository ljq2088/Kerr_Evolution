"""Isolated Li-prefactor order-four Green function at real positive frequency."""
import numpy as np
from scipy.integrate import solve_ivp
from environment_cloud import angular_eigenvalue,radial_coefficients
from environment_radial import horizon_series
from report_li_order4_boundary_control import formal_boundary

class LiOrder4Green:
    def __init__(self,a,mu,omega,ell,m,rmax=None,offset=1e-4,rtol=1e-11):
        if omega<=0:raise ValueError("This driver is restricted to positive-frequency channels")
        if abs(omega**2-mu**2)<1e-10:raise ValueError("Unresolved threshold")
        self.a,self.mu,self.omega,self.m=a,mu,omega,m
        self.rp=1+np.sqrt(1-a*a);self.rmin=self.rp+offset
        self.propagating=bool(omega>mu)
        self.rmax=float(rmax if rmax is not None else (8000. if self.propagating else 1000.))
        if self.rmax<=self.rmin or offset<=0:raise ValueError("Invalid domain")
        self.lam=angular_eigenvalue(ell,m,a*a*(omega*omega-mu*mu))
        def rhs(r,y):
            d,dp,v=radial_coefficients(r,a,mu,omega,m,self.lam)
            return [y[1],-(dp*y[1]+v*y[0])/d]
        start=horizon_series(self.rmin,a,mu,omega,m,self.lam,order=4)
        up,dup,audit=formal_boundary(a,mu,omega,m,self.lam,self.rmax,order=4,dps=64,propagating=self.propagating)
        self.boundary_audit=audit
        opts=dict(method='DOP853',dense_output=True,rtol=rtol,atol=rtol*1e-3)
        self.insol=solve_ivp(rhs,(self.rmin,self.rmax),start,**opts)
        self.upsol=solve_ivp(rhs,(self.rmax,self.rmin),[up,dup],**opts)
        if not self.insol.success or not self.upsol.success:raise RuntimeError("Radial integration failed")
        self.w0=self.wronskian(np.array([np.sqrt(self.rmin*self.rmax)]))[0]
        if not np.isfinite(self.w0) or self.w0==0:raise ValueError("Invalid Wronskian")
    def wronskian(self,r):
        u,du=self.insol.sol(r);v,dv=self.upsol.sol(r)
        return (r*r-2*r+self.a*self.a)*(u*dv-v*du)
