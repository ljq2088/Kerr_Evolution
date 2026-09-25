"""Leaver cloud with explicit finite-BL mass and temporal-only freezing.

The author's horizon prescription is not public. This adapter records rather
than hides its finite BL t=0 integration limits; it does not assert equality
with an integral extending to the mathematical horizon.
"""
import numpy as np
from scipy.integrate import simpson
from li_leaver_cloud import LiLeaverCloud
from environment_kerr_quasibound_cloud import ComplexScalarAngular

class LiNormalizedCloud:
    def __init__(self,ell=1,terms=150,dps=64,mass_offset=1e-6,mass_outer=2000.,mass_points=12001):
        if mass_offset<=0 or mass_outer<=2 or mass_points<101:
            raise ValueError("Invalid finite mass domain")
        self.leaver=LiLeaverCloud(ell=ell,m=ell,terms=terms,dps=dps)
        c=self.leaver
        self.a,self.mu,self.m,self.ell=.88,.3,ell,ell
        self.spectral_omega=complex(c.omega)
        self.omega=self.spectral_omega.real
        self.rp,self.rm=float(c.data['rp']),float(c.data['rm'])
        self.rmax=mass_outer
        self.angular_eigenfunction=ComplexScalarAngular(ell,ell,self.a**2*(self.spectral_omega**2-self.mu**2))
        self.lam=self.angular_eigenfunction.A
        self.coeff=np.array([complex(v) for v in c.coefficients])
        self.q,self.beta,self.s=[complex(c.data[key]) for key in ('q','beta','s')]
        self.amplitude=1.
        raw=self.bl_integrals(mass_offset,mass_outer,mass_points)
        if not np.isfinite(raw['energy']) or raw['energy']<=0:
            raise ValueError("Finite BL mass is not positive")
        self.amplitude=1/np.sqrt(raw['energy'])
        self.provenance=dict(a=self.a,mu=self.mu,ell_c=ell,m_c=ell,n_c=0,
            spectral_omega=[self.spectral_omega.real,self.spectral_omega.imag],
            temporal_omega=self.omega,frequency_policy='complex spatial Leaver profile; temporal frequency Re(omega_c)',
            radial_terms=terms,spectral_decimal_digits=dps,evaluation_precision='complex128',
            mass_prescription='finite BL t=0 Killing energy using frozen temporal frequency',
            mass_inner_offset=mass_offset,mass_outer=mass_outer,mass_points=mass_points,
            raw_bl_mass=raw['energy'],amplitude=self.amplitude,
            normalized_bl_mass=1.,author_horizon_prescription_verified=False)

    def radial_state(self,r):
        r=np.atleast_1d(np.asarray(r,float))
        if np.any(r<=self.rp):raise ValueError("Exterior radii required")
        gap=self.rp-self.rm;x=(r-self.rp)/(r-self.rm)
        xp=gap/(r-self.rm)**2;xpp=-2*gap/(r-self.rm)**3
        p=np.polynomial.polynomial
        P=p.polyval(x,self.coeff);Px=p.polyval(x,p.polyder(self.coeff))
        Pxx=p.polyval(x,p.polyder(self.coeff,2))
        F=self.amplitude*np.exp(self.q*(r-self.rp)+self.beta*np.log((r-self.rm)/gap)+self.s*np.log(x))
        L=self.q+self.beta/(r-self.rm)+self.s*xp/x
        Lp=-self.beta/(r-self.rm)**2+self.s*(xpp/x-(xp/x)**2)
        return np.array([F*P,F*(L*P+Px*xp),F*((L*L+Lp)*P+2*L*Px*xp+Pxx*xp*xp+Px*xpp)])

    def radial(self,r):return self.radial_state(r)[:2]

    def angular(self,theta):
        s,sp=self.angular_eigenfunction.evaluate(theta)
        return s,sp,self.lam

    def angular_state(self,theta):
        s,sp,A=self.angular(theta);st=np.sin(theta);ct=np.cos(theta)
        c2=self.a**2*(self.spectral_omega**2-self.mu**2)
        return np.array([s,sp,-ct/st*sp-(c2*ct**2-self.m**2/st**2+A)*s])

    def bl_integrals(self,offset,outer,nr=12001):
        r=self.rp+np.geomspace(offset,outer-self.rp,nr)
        R,Rp=self.radial(r)
        x,weights=np.polynomial.legendre.leggauss(48)
        S,Sp,_=self.angular(np.arccos(x))
        rr=r[:,None];z=x[None,:]
        sig=rr*rr+self.a**2*z*z;d=(rr-self.rp)*(rr-self.rm)
        gtt=-((rr*rr+self.a**2)**2-self.a**2*d*(1-z*z))/(sig*d)
        gtp=-2*self.a*rr/(sig*d)
        gpp=(d-self.a**2*(1-z*z))/(sig*d*(1-z*z))
        f2=abs(R[:,None])**2*abs(S)**2
        E=sig*((-gtt*self.omega**2+gpp*self.m**2+self.mu**2)*f2)
        E+=d*abs(Rp[:,None])**2*abs(S)**2+abs(R[:,None])**2*abs(Sp)**2
        Q=2*sig*(-gtt*self.omega+gtp*self.m)*f2
        return dict(energy=float(2*np.pi*simpson(E@weights,x=r)),
                    charge=float(2*np.pi*simpson(Q@weights,x=r)))

    def normalization_audit(self):
        return [dict(inner_offset=e,outer_radius=o,**self.bl_integrals(e,o))
                for e,o in [(1e-3,2000.),(1e-6,1000.),(1e-6,2000.),(1e-9,2000.)]]

    def hessian(self,r,theta):
        from environment_source import connection
        R,Rp,Rpp=self.radial_state(r)[:,0]
        S,Sp,Spp=self.angular_state(np.asarray(theta))
        field=R*S
        grad=np.array([-1j*self.omega*field,Rp*S,R*Sp,1j*self.m*field])
        partial=np.empty((4,4),complex)
        partial[0]=-1j*self.omega*grad;partial[:,0]=partial[0]
        partial[3]=1j*self.m*grad;partial[:,3]=partial[3]
        partial[1,1],partial[2,2]=Rpp*S,R*Spp
        partial[1,2]=partial[2,1]=Rp*Sp
        inverse,gamma=connection(r,theta,self.a)
        return field,partial-np.einsum('kij,k->ij',gamma,grad),inverse

    def lorenz_source(self,r,theta,h_covariant):
        _,hessian,inverse=self.hessian(r,theta)
        return np.einsum('ij,ij->',inverse@h_covariant@inverse,hessian)
