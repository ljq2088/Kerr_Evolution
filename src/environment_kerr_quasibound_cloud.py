"""Bounded general-spin |211> complex cloud, normalized on ingoing KS T=0.

This new module does not replace production's synchronized real cloud. A
supplied complex spectrum must be an actual quasibound root. The default root
is independently Leaver/shooting-validated at a=.88, mu=.3. No spectral solver
is silently applied to other parameters.

freeze_growth=True is an explicit temporal-only approximation AFTER the exact
complex spatial eigenfunction and its KS-slice mass normalization are formed.
"""
import numpy as np
from scipy.integrate import solve_ivp, simpson
from scipy.special import lpmv,gammaln
from environment_cloud import radial_coefficients
from environment_radial import infinity_series
from environment_source import connection

FIXED_A088_OMEGA=complex(.296293534725611462814885674586,2.2166093964457949473e-9)

class ComplexScalarAngular:
    """Complex-symmetric angular eigenproblem; physical Hermitian norm=1.

    Distinct complex-frequency eigenfunctions are orthogonal bilinearly, not
    Hermitian-orthogonal. dual(theta)=S(theta)/sum(coeff**2) is the projection
    weight in theta (the azimuthal Fourier coefficient is already removed).
    """
    def __init__(self,ell,m,c2,size=24):
        if int(ell)!=ell or int(m)!=m or ell<abs(m):
            raise ValueError("Integer ell >= |m| required")
        self.ell,self.m,self.c2=int(ell),int(m),complex(c2)
        mm=abs(self.m);self.ls=np.arange(mm,mm+max(size,ell-mm+12))
        ext=np.r_[self.ls,self.ls[-1]+1];cosine=np.zeros((len(ext),len(ext)))
        for j,l in enumerate(ext[:-1]):
            cosine[j,j+1]=cosine[j+1,j]=np.sqrt(((l+1)**2-mm*mm)/((2*l+1)*(2*l+3)))
        matrix=np.diag(self.ls*(self.ls+1.))-self.c2*(cosine@cosine)[:-1,:-1]
        values,vectors=np.linalg.eig(matrix)
        index=int(np.argmin(abs(values-ell*(ell+1))))
        coeff=vectors[:,index];pivot=ell-mm
        coeff=coeff*np.exp(-1j*np.angle(coeff[pivot]))/np.linalg.norm(coeff)
        self.A=complex(values[index]);self.coefficients=coeff
        self.bilinear_norm=complex(np.dot(coeff,coeff))
        self.eigenpair_residual=float(np.linalg.norm(matrix@coeff-self.A*coeff))
        if abs(self.bilinear_norm)<1e-12:
            raise ValueError("Near self-orthogonal angular mode; dual projection ill-conditioned")
    def evaluate(self,theta):
        theta=np.asarray(theta,float);x=np.cos(theta);st=np.sin(theta);mm=abs(self.m)
        if np.any(theta<=0) or np.any(theta>=np.pi) or not np.all(np.isfinite(theta)):
            raise ValueError("Use finite polar angles strictly between the axes")
        value=np.zeros(theta.shape,complex);slope=np.zeros(theta.shape,complex)
        sign=(-1)**mm if self.m<0 else 1
        for l,b in zip(self.ls,self.coefficients):
            norm=np.sqrt((2*l+1)/(4*np.pi)*np.exp(gammaln(l-mm+1)-gammaln(l+mm+1)))
            P=lpmv(mm,l,x);Pm=lpmv(mm,l-1,x) if l>mm else np.zeros_like(x)
            value+=sign*b*norm*P
            slope+=sign*b*norm*(l*x*P-(l+mm)*Pm)/st
        return value,slope
    def dual(self,theta):
        return self.evaluate(theta)[0]/self.bilinear_norm

class GeneralKerrQuasiboundCloud:
    """BL source interface backed by a horizon-regular KS radial integration."""
    def __init__(self,a=.88,alpha=.3,omega=None,mass=1.,ell=1,m=1,
                 freeze_growth=False,offset=1e-6,outer_efolds=70.,
                 rtol=2e-12,radial_points=6001,angular_size=24,
                 horizon_order=6,infinity_order=8,match_radius=None):
        if not 0<=a<1 or alpha<=0 or mass<=0 or offset<=0:
            raise ValueError("Require 0<=a<1 and positive mu, mass, horizon offset")
        if (ell,m)!=(1,1):
            raise ValueError("This bounded implementation is validated for |211> only")
        if omega is None:
            if a!=.88 or alpha!=.3:
                raise ValueError("Supply a independently established complex quasibound omega")
            omega=FIXED_A088_OMEGA
        self.a,self.mu,self.ell,self.m=float(a),float(alpha),ell,m
        self.spectral_omega=complex(omega);self.freeze_growth=bool(freeze_growth)
        if not np.isfinite(self.spectral_omega) or not 0<self.spectral_omega.real<alpha:
            raise ValueError("Require finite positive-frequency quasibound spectrum with Re(omega)<mu")
        self.omega=self.spectral_omega.real if freeze_growth else self.spectral_omega
        self.rp=1+np.sqrt(1-a*a);self.rm=2-self.rp
        self.c2=a*a*(self.spectral_omega**2-alpha**2)
        self.angular_eigenfunction=ComplexScalarAngular(ell,m,self.c2,angular_size)
        self.lam=self.angular_eigenfunction.A
        self.rmin=self.rp+offset
        self.decay_kappa=np.sqrt(alpha*alpha-self.spectral_omega**2)
        if self.decay_kappa.real<=0:raise ValueError("Need decaying infinity branch")
        self.rmax=float(outer_efolds/self.decay_kappa.real)
        self.match=2/alpha**2 if match_radius is None else float(match_radius)
        if not self.rmin<self.match<self.rmax:raise ValueError("Invalid matching interval")
        w=self.spectral_omega;d=self.rp-self.rm
        self.horizon_KS_unit=np.exp(-1j*w*self.rp+1j*m*(a/2+a/self.rp*np.log(d/2)))
        def rhs_ks(r,state):
            delta=(r-self.rp)*(r-self.rm)
            B=2*(r-1)-2j*(2*r*w-a*m)
            V=(w*w-alpha*alpha)*r*r+2*r*w*w-self.lam-2j*w
            return [state[1],-(B*state[1]+V*state[0])/delta]
        self.rhs_ks=rhs_ks
        def rhs_bl(r,state):
            delta,dp,V=radial_coefficients(r,a,alpha,w,m,self.lam)
            return [state[1],-(dp*state[1]+V*state[0])/delta]
        self.rhs=rhs_bl
        b0=d-2j*(2*self.rp*w-a*m);b1=2-4j*w
        v0=(w*w-alpha*alpha)*self.rp**2+2*self.rp*w*w-self.lam-2j*w
        v1=2*self.rp*(w*w-alpha*alpha)+2*w*w;v2=w*w-alpha*alpha
        coefficients=[complex(self.horizon_KS_unit)]
        for n in range(1,horizon_order+1):
            known=((n-1)*(n-2)+b1*(n-1)+v0)*coefficients[n-1]
            if n>=2:known+=v1*coefficients[n-2]
            if n>=3:known+=v2*coefficients[n-3]
            coefficients.append(-known/(n*((n-1)*d+b0)))
        self.ks_horizon_series=np.asarray(coefficients)
        F=sum(c*offset**n for n,c in enumerate(coefficients))
        Fp=sum(n*c*offset**(n-1) for n,c in enumerate(coefficients) if n)
        k=1j*self.decay_kappa
        power,f,df=infinity_series(self.rmax,a,alpha,w,m,self.lam,k,order=infinity_order)
        phase,H=self._phase(self.rmax,w)
        options=dict(method="DOP853",rtol=rtol,atol=rtol*1e-3,dense_output=True)
        self.left=solve_ivp(rhs_ks,(self.rmin,self.match),[F,Fp],**options)
        self.right=solve_ivp(rhs_ks,(self.rmax,self.match),
            [phase,phase*(1j*k+power/self.rmax+df/f+1j*H)],**options)
        if not self.left.success or not self.right.success:raise RuntimeError("KS radial integration failed")
        left,right=self.left.y[:,-1],self.right.y[:,-1]
        self.glue=left[0]/right[0]
        self.matching_log_derivative_residual=complex(left[1]/left[0]-right[1]/right[0])
        if abs(self.matching_log_derivative_residual)>1e-8:
            raise ValueError("Supplied spectrum/boundaries fail quasibound derivative matching")
        self.amplitude=1.
        E,Q=self.ks_integrals(radial_points)
        if not np.isfinite(E) or E<=0:raise ValueError("Nonpositive or nonfinite Killing energy")
        self.amplitude=float(np.sqrt(mass/E))
        self.mass,self.charge=mass,Q*self.amplitude**2
        self.normalization="Underlying exact complex mode Killing energy on ingoing Kerr-Schild T=0"
        self.normalization_frequency=self.spectral_omega
        self.normalization_audit=dict(raw_energy=E,raw_charge=Q,radial_points=radial_points,
            angular_order=48,offset=offset,outer_efolds=outer_efolds,rtol=rtol,
            horizon_order=horizon_order,infinity_order=infinity_order,
            freeze_order="Time frequency only, AFTER complex spatial eigenfunction and mass normalization" if freeze_growth else "No frequency freezing")
    def _phase(self,r,w):
        r=np.asarray(r);d=self.rp-self.rm
        h=2*self.rp/d*np.log((r-self.rp)/2)-2*self.rm/d*np.log((r-self.rm)/2)
        g=self.a/d*np.log((r-self.rp)/(r-self.rm))
        delta=(r-self.rp)*(r-self.rm)
        return np.exp(1j*w*h-1j*self.m*g),(2*r*w-self.a*self.m)/delta
    def _check_r(self,r):
        r=np.atleast_1d(r).astype(float)
        if not np.all(np.isfinite(r)) or np.any(r<self.rmin) or np.any(r>self.rmax):
            raise ValueError(f"Requested cloud radius outside solved [{self.rmin}, {self.rmax}]")
        return r
    def ks_radial(self,r,frequency=None):
        """Underlying spectral mode F,F'; optional different temporal frequency is explicit."""
        r=self._check_r(r);state=np.empty((2,len(r)),complex);mask=r<=self.match
        if mask.any():state[:,mask]=self.left.sol(r[mask])
        if (~mask).any():state[:,~mask]=self.glue*self.right.sol(r[~mask])
        state*=self.amplitude
        if frequency is not None and complex(frequency)!=self.spectral_omega:
            old,Ho=self._phase(r,self.spectral_omega);new,Hn=self._phase(r,frequency)
            state=np.array([new/old*state[0],new/old*(state[1]+1j*(Hn-Ho)*state[0])])
        return state
    def radial(self,r):
        r=self._check_r(r);F,Fp=self.ks_radial(r);phase,H=self._phase(r,self.spectral_omega)
        return np.array([F/phase,(Fp-1j*H*F)/phase])
    def radial_state(self,r):
        r=self._check_r(r);R,Rp=self.radial(r)
        return np.array([R,Rp,self.rhs(r,[R,Rp])[1]])
    def angular(self,theta):
        S,Sp=self.angular_eigenfunction.evaluate(theta)
        return S,Sp,self.lam
    def hessian(self,r,theta):
        R,Rp,Rpp=self.radial_state(r)[:,0];S,Sp=self.angular_eigenfunction.evaluate(theta)
        Spp=-np.cos(theta)/np.sin(theta)*Sp-(self.c2*np.cos(theta)**2-self.m**2/np.sin(theta)**2+self.lam)*S
        field=R*S
        grad=np.array([-1j*self.omega*field,Rp*S,R*Sp,1j*self.m*field])
        partial=np.empty((4,4),complex)
        partial[0]=-1j*self.omega*grad;partial[:,0]=partial[0]
        partial[3]=1j*self.m*grad;partial[:,3]=partial[3]
        partial[1,1],partial[2,2]=Rpp*S,R*Spp
        partial[1,2]=partial[2,1]=Rp*Sp
        inverse,gamma=connection(r,theta,self.a)
        return field,partial-np.einsum("kij,k->ij",gamma,grad),inverse
    def lorenz_source(self,r,theta,h_covariant):
        _,hessian,inverse=self.hessian(r,theta);h=np.asarray(h_covariant,complex)
        if h.shape!=(4,4) or not np.allclose(h,h.T):raise ValueError("Need symmetric covariant BL metric")
        return np.einsum("ij,ij->",inverse@h@inverse,hessian)
    def jet(self,r,theta,order=6,return_geometry=False):
        from lorenz_ghp import KerrGHP
        from lorenz_mode_jet import separated_jet
        spectral=KerrGHP(r,theta,self.a,omega=self.spectral_omega,m=self.m,order=order)
        R,Rp=self.radial(r)[:,0];S,Sp=self.angular_eigenfunction.evaluate(theta)
        radial_lambda=self.lam+self.a**2*self.spectral_omega**2-2*self.a*self.m*self.spectral_omega
        field=separated_jet(spectral,0,radial_lambda,R,Rp,S,Sp,mass_squared=self.mu**2)
        geometry=spectral
        if self.freeze_growth:
            geometry=KerrGHP(r,theta,self.a,omega=self.omega,m=self.m,order=order)
        return (geometry,field) if return_geometry else field
    def ks_integrals(self,radial_points=6001,ntheta=48,frequency=None):
        """Exact stress-energy/Noether integrals, not the shortcut E=omega Q.

        Default always normalizes the underlying spectral solution. Supplying
        frequency=self.omega diagnoses the temporally frozen approximate field.
        """
        r=self.rp+np.geomspace(self.rmin-self.rp,self.rmax-self.rp,radial_points)
        r[0],r[-1]=self.rmin,self.rmax
        w=self.spectral_omega if frequency is None else complex(frequency)
        F,Fp=self.ks_radial(r,frequency=w)
        x,wt=np.polynomial.legendre.leggauss(ntheta);S,Sp=self.angular_eigenfunction.evaluate(np.arccos(x))
        rr=r[:,None];sig=rr**2+self.a**2*x**2;delta=(rr-self.rp)*(rr-self.rm)
        fn=abs(F[:,None])**2;sn=abs(S)**2;cross=np.imag(F.conjugate()*Fp)[:,None]
        energy=((sig+2*rr)*abs(w)**2+self.mu**2*sig+self.m**2/(1-x*x))*fn*sn
        energy+=delta*abs(Fp[:,None])**2*sn+fn*abs(Sp)**2+2*self.a*self.m*cross*sn
        charge=2*((sig+2*rr)*w.real*fn+2*rr*cross)*sn
        return float(2*np.pi*simpson(energy@wt,x=r)),float(2*np.pi*simpson(charge@wt,x=r))
    def integrals(self,radial_points=6001):
        return self.ks_integrals(radial_points)
    @property
    def provenance(self):
        import hashlib
        from pathlib import Path
        return dict(kind="GeneralKerrQuasiboundCloud",a=self.a,mu=self.mu,ell=self.ell,m=self.m,
            spectral_omega=[self.spectral_omega.real,self.spectral_omega.imag],
            temporal_omega=[complex(self.omega).real,complex(self.omega).imag],
            freeze_growth=self.freeze_growth,mass=self.mass,normalization=self.normalization,
            normalization_audit=self.normalization_audit,
            radial_coordinates="BL",normalization_slice="ingoing Kerr-Schild T=0",
            azimuth_transform="Psi=phi+a/(rplus-rminus)*log((r-rplus)/(r-rminus))",
            phase_convention="Unit exp[-i(omega-m OmegaH) rstar] ingoing BL amplitude before mass normalization",
            module_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    def horizon_balance(self,radial_points=6001):
        E,Q=self.ks_integrals(radial_points);w=self.spectral_omega
        norm=abs(self.amplitude*self.horizon_KS_unit)**2
        FE=2*(2*self.rp*abs(w)**2-self.a*self.m*w.real)*norm
        FQ=2*(2*self.rp*w.real-self.a*self.m)*norm
        return dict(energy=E,charge=Q,horizon_energy=float(FE),horizon_charge=float(FQ),
            energy_balance=float(2*w.imag*E+FE),charge_balance=float(2*w.imag*Q+FQ),
            energy_balance_relative=float(abs(2*w.imag*E+FE)/max(abs(FE),1e-300)),
            charge_balance_relative=float(abs(2*w.imag*Q+FQ)/max(abs(FQ),1e-300)))

def project_complex_source(cloud,radii,orbital_radius,ell,m,metric_mode,ntheta=48):
    """Complex-frequency separated source using the bilinear angular dual.

    For real frequency this agrees with the old real angular projection.
    Does not modify production project_source or claim real-frequency flux
    formulas apply to a growing complex-frequency forced mode.
    """
    omega=cloud.omega+(m-cloud.m)/(orbital_radius**1.5+cloud.a)
    angular=ComplexScalarAngular(ell,m,cloud.a**2*(omega**2-cloud.mu**2))
    x,weights=np.polynomial.legendre.leggauss(ntheta);theta=np.arccos(x)
    dual=angular.dual(theta);values=[]
    for r in radii:
        source=np.array([cloud.lorenz_source(r,t,metric_mode(r,t)) for t in theta])
        values.append(2*np.pi*np.dot(weights,dual*(r*r+cloud.a**2*x*x)*source))
    return omega,np.asarray(values)
