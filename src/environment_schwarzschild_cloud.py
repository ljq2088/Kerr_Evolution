"""Schwarzschild massive scalar |211> quasibound spectrum, M=1.

The complex eigenfrequency is retained here. The paper's later replacement
Im(omega)=0 in the environmental source is a separate approximation.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from environment_cloud import radial_coefficients
from environment_radial import horizon_series,infinity_series
from environment_source import ThresholdCloud
from scipy.integrate import simpson


def quasibound_211(mu=.3,outer_efolds=45.,offset=1e-5,rtol=2e-11,seed=None):
    if not .19<=mu<=.31:
        raise ValueError('Only the paper alpha=.2/.3 neighborhood is targeted')
    rmin=2+offset
    rmax=2*outer_efolds/mu**2
    match=2/mu**2
    if not 2<rmin<match<rmax:
        raise ValueError('Invalid shooting intervals')
    def integrate(omega):
        k=1j*np.sqrt(mu*mu-omega*omega+0j)
        if k.imag<=0:
            raise ValueError('Must select the infinity-decaying branch')
        def rhs(r,y):
            d,dp,v=radial_coefficients(r,0.,mu,omega,1,2.)
            return [y[1],-(dp*y[1]+v*y[0])/d]
        rin,din=horizon_series(rmin,0.,mu,omega,1,2.)
        p,f,df=infinity_series(rmax,0.,mu,omega,1,2.,k)
        left=solve_ivp(rhs,(rmin,match),[rin,din],method='DOP853',rtol=rtol,atol=rtol*1e-3)
        right=solve_ivp(rhs,(rmax,match),[1.+0j,1j*k+p/rmax+df/f],method='DOP853',rtol=rtol,atol=rtol*1e-3)
        if not left.success or not right.success:
            raise RuntimeError('Quasibound radial shooting failed')
        u,v=left.y[:,-1],right.y[:,-1]
        residual=u[1]/u[0]-v[1]/v[0]
        return residual
    if seed is None:
        seed=mu*(1-mu*mu/8)-1j*(1e-6 if mu>.25 else 1e-8)
    def equations(x):
        value=integrate(complex(*x))
        return [value.real,value.imag]
    result=root(equations,[seed.real,seed.imag],tol=1e-12)
    omega=complex(*result.x)
    # The damping rate is much smaller than Re(omega); a relative stopping
    # test on the two-component root alone can leave Im(omega) inaccurate.
    for polish in range(6):
        residual=integrate(omega)
        if abs(residual)<1e-13:
            break
        step=1e-6
        derivative=(integrate(omega+step)-integrate(omega-step))/(2*step)
        omega-=residual/derivative
    residual=integrate(omega)
    if abs(residual)>1e-11 or not 0<omega.real<mu or omega.imag>=0:
        raise RuntimeError(f'Invalid quasibound root: {result.message}; omega={omega}; residual={residual}')
    return dict(alpha=mu,ell=1,m=1,a_over_M=0.,
        omega=[omega.real,omega.imag],log_derivative_residual=[residual.real,residual.imag],
        outer_efolds=outer_efolds,horizon_offset=offset,rmax=rmax,match_radius=match,
        solver_success=bool(result.success),function_evaluations=int(result.nfev),
        newton_polish_iterations=polish,
        amplitude_decay_time_M=-1/omega.imag)


class SchwarzschildCloud(ThresholdCloud):
    """Quasibound BL profile normalized on ingoing Kerr-Schild T=0.

    T=t+2 log((r-2)/2). Setting freeze_decay=True changes temporal frequency
    only, explicitly reproducing a non-exact stationary-background prescription.
    It does not change the underlying complex-frequency radial profile.
    """
    def __init__(self,alpha=.3,mass=1.,freeze_decay=False,offset=1e-6,radial_points=6001):
        if mass<=0:
            raise ValueError('Cloud mass must be positive')
        spectrum=quasibound_211(alpha,outer_efolds=60.,offset=offset,rtol=2e-12)
        self.spectral_omega=complex(*spectrum['omega'])
        self.omega=self.spectral_omega.real if freeze_decay else self.spectral_omega
        self.a,self.mu,self.m,self.rp,self.c2,self.lam=0.,alpha,1,2.,0.,2.
        self.rmin,self.rmax,self.match=2+offset,spectrum['rmax'],spectrum['match_radius']
        self.freeze_decay=freeze_decay
        self.normalization='Killing energy of complex mode on ingoing Kerr-Schild T=0'
        w=self.spectral_omega
        def rhs(r,y):
            d,dp,v=radial_coefficients(r,0.,alpha,w,1,2.)
            return [y[1],-(dp*y[1]+v*y[0])/d]
        self.rhs=rhs
        rin,din=horizon_series(self.rmin,0.,alpha,w,1,2.)
        k=1j*np.sqrt(alpha*alpha-w*w)
        p,f,df=infinity_series(self.rmax,0.,alpha,w,1,2.,k)
        options=dict(method='DOP853',rtol=2e-12,atol=2e-15,dense_output=True)
        self.left=solve_ivp(rhs,(self.rmin,self.match),[rin,din],**options)
        self.right=solve_ivp(rhs,(self.rmax,self.match),[1.+0j,1j*k+p/self.rmax+df/f],**options)
        if not self.left.success or not self.right.success:
            raise RuntimeError('Schwarzschild cloud profile integration failed')
        self.glue=self.left.y[0,-1]/self.right.y[0,-1]
        self.amplitude=1.
        energy,charge=self.ks_integrals(radial_points)
        self.amplitude=np.sqrt(mass/energy)
        self.mass,self.charge=mass,charge*self.amplitude**2

    def radial(self,r):
        r=np.atleast_1d(r).astype(float)
        if np.any(r<self.rmin) or np.any(r>self.rmax):
            raise ValueError('Outside cloud radial domain')
        values=np.empty((2,r.size),complex)
        left=r<=self.match
        if left.any():values[:,left]=self.left.sol(r[left])
        if (~left).any():values[:,~left]=self.glue*self.right.sol(r[~left])
        return self.amplitude*values

    def ks_radial(self,r):
        r=np.atleast_1d(r).astype(float)
        R,dR=self.radial(r)
        phase=np.exp(2j*self.spectral_omega*np.log((r-2)/2))
        return phase*R,phase*(dR+2j*self.spectral_omega*R/(r-2))

    def ks_integrals(self,radial_points=6001):
        r=2+np.geomspace(self.rmin-2,self.rmax-2,radial_points)
        r[0],r[-1]=self.rmin,self.rmax
        R,dR=self.ks_radial(r)
        w=self.spectral_omega
        energy=r*r*((1+2/r)*abs(w)**2*abs(R)**2+(1-2/r)*abs(dR)**2
                    +(2/r**2+self.mu**2)*abs(R)**2)
        charge=2*r*r*((1+2/r)*w.real*abs(R)**2+2/r*np.imag(R.conjugate()*dR))
        return float(simpson(energy,x=r)),float(simpson(charge,x=r))

    def decay_balance(self):
        E,Q=self.ks_integrals()
        Rh=self.ks_radial(self.rmin)[0][0]
        horizon_energy=8*abs(self.spectral_omega)**2*abs(Rh)**2
        horizon_charge=8*self.spectral_omega.real*abs(Rh)**2
        return dict(energy=E,charge=Q,horizon_energy=float(horizon_energy),
                    horizon_charge=float(horizon_charge),
                    energy_balance=float(2*self.spectral_omega.imag*E+horizon_energy),
                    charge_balance=float(2*self.spectral_omega.imag*Q+horizon_charge))

    def integrals(self,radial_points=6001):
        """Use the horizon-penetrating slice, not the inherited BL integral."""
        return self.ks_integrals(radial_points)


def main():
    rows=[]
    for alpha in (.2,.3):
        baseline=quasibound_211(alpha)
        refined=quasibound_211(alpha,outer_efolds=60.,offset=1e-6,rtol=2e-12,
                              seed=complex(*baseline['omega']))
        delta=abs(complex(*baseline['omega'])-complex(*refined['omega']))
        rows.append(dict(baseline=baseline,refined=refined,frequency_cutoff_change=delta))
        print(json.dumps(rows[-1]),flush=True)
    output=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/schwarzschild_cloud_spectrum.json'
    output.write_text(json.dumps(dict(status='complex_background_spectrum_only_not_schwarzschild_flux_reproduction',
        method='Independent two-sided radial shooting with ingoing/decaying boundaries',runs=rows),indent=2)+'\n')
    balances=[]
    for alpha in (.2,.3):
        cloud=SchwarzschildCloud(alpha=alpha)
        balances.append(dict(alpha=alpha,omega=[cloud.spectral_omega.real,cloud.spectral_omega.imag],
                             balance=cloud.decay_balance()))
    output.with_name('schwarzschild_cloud_balance.json').write_text(json.dumps(dict(
        status='exact_complex_cloud_KS_slice_validation',runs=balances),indent=2)+'\n')


if __name__=='__main__':
    main()
