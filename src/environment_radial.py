"""Massive Kerr radial In/Up Green solver, M=1, real nonthreshold frequency.

Boundary expansions retain 4 horizon and 6 infinity corrections. Cutoff convergence is required;
near k=0 and narrow resonances a dedicated asymptotic treatment is still needed.
"""
import numpy as np
from scipy.integrate import solve_ivp, cumulative_simpson
from environment_cloud import angular_eigenvalue, radial_coefficients


def tortoise(r, a):
    rp, rm = 1+np.sqrt(1-a*a), 1-np.sqrt(1-a*a)
    return r+2*rp/(rp-rm)*np.log((r-rp)/2)-2*rm/(rp-rm)*np.log((r-rm)/2)


def complex_cumulative(y, x):
    return (cumulative_simpson(y.real, x=x, initial=0)+
            1j*cumulative_simpson(y.imag, x=x, initial=0))


def infinity_series(r, a, mu, omega, m, lam, k, order=6):
    """Outgoing/decaying e^(ik r) r^p sum c_n/r^n; solve polynomial recurrence."""
    p = -1+1j*(2*omega*omega-mu*mu)/k
    d2 = {4:1.,3:-4.,2:4+2*a*a,1:-4*a*a,0:a**4}
    dd = {3:2.,2:-6.,1:4+2*a*a,0:-2*a*a}
    C, b = a*a*omega*omega-2*a*m*omega+lam, a*a*omega-a*m
    potential = {4:k*k,3:2*mu*mu,2:2*omega*b-C-a*a*mu*mu,1:2*C,0:b*b-a*a*C}
    def coefficient(n, target):
        total = 0j
        for degree,v in d2.items():
            for shift,f in ((0,-k*k),(1,2j*k*(p-n)),(2,(p-n)*(p-n-1))):
                if degree-n-shift == target:
                    total += v*f
        for degree,v in dd.items():
            for shift,f in ((0,1j*k),(1,p-n)):
                if degree-n-shift == target:
                    total += v*f
        total += potential.get(target+n,0.)
        return total
    coeff = [1.+0j]
    for n in range(1,order+1):
        target=3-n
        coeff.append(-sum(c*coefficient(j,target) for j,c in enumerate(coeff))/coefficient(n,target))
    f=sum(c/r**n for n,c in enumerate(coeff))
    df=sum(-n*c/r**(n+1) for n,c in enumerate(coeff))
    return p,f,df


def horizon_series(r, a, mu, omega, m, lam, order=4):
    """Ingoing Frobenius series with unit exp[-i(omega-m OmegaH)r*]."""
    from numpy.polynomial import Polynomial as P
    rp,rm = 1+np.sqrt(1-a*a),1-np.sqrt(1-a*a)
    d,x = rp-rm,r-rp
    rr=P([rp,1.])
    delta=P([0.,d,1.])
    d2,dd=delta*delta,delta*delta.deriv()
    potential=((rr*rr+a*a)*omega-a*m)**2-delta*(mu*mu*rr*rr+a*a*omega*omega-2*a*m*omega+lam)
    dw=omega-m*a/(2*rp)
    s=-1j*dw*2*rp/d
    def coefficient(n,target):
        value=0j
        for j,v in enumerate(d2.coef):
            if n-2+j==target:
                value+=v*(s+n)*(s+n-1)
        for j,v in enumerate(dd.coef):
            if n-1+j==target:
                value+=v*(s+n)
        for j,v in enumerate(potential.coef):
            if n+j==target:
                value+=v
        return value
    coeff=[1.+0j]
    for n in range(1,order+1):
        coeff.append(-sum(c*coefficient(j,n) for j,c in enumerate(coeff))/coefficient(n,n))
    f=sum(c*x**n for n,c in enumerate(coeff))
    df=sum(n*c*x**(n-1) for n,c in enumerate(coeff) if n)
    constant=rp-2*rp/d*np.log(2)-2*rm/d*np.log(d/2)
    value=np.exp(s*np.log(x)-1j*dw*constant)*f
    return value,value*(s/x+df/f)


def coulomb_boundary(r,a,mu,omega,lam,k,propagating):
    """Whittaker boundary for y=sqrt(Delta) R, retaining Q through r^-2.

    y''+[k²+2 beta/r+C2/r²+O(r^-3)]y=0, beta=2 omega²-mu².
    Unlike the inverse-r wave series this resums beta/(k² r) near threshold.
    The neglected O(r^-3) potential still requires outer-cutoff convergence.
    """
    import mpmath as mp
    with mp.workdps(40):
        kk=mp.mpc(k)
        beta=2*omega*omega-mu*mu
        c2=12*omega*omega-4*mu*mu-a*a*k*k-lam
        kap=1j*mp.mpc(beta)/kk
        nu=mp.sqrt(mp.mpc(.25-c2))
        z=-2j*kk*r
        value=mp.whitw(kap,nu,z)
        derivative=mp.diff(lambda zz:mp.whitw(kap,nu,zz),z)*(-2j*kk)
        delta=r*r-2*r+a*a
        logarithmic=derivative/value-(r-1)/delta
        if propagating:
            norm=mp.exp(-kap*mp.log(-2j*kk)-2j*kk*mp.log(2))
            up=complex(norm*value/mp.sqrt(delta))
        else:
            up=1.+0j
        return up,up*complex(logarithmic)


class RadialGreen:
    def __init__(self, a, mu, omega, ell, m, rmax=2000., offset=1e-6, rtol=2e-10,
                 mass_squared=None,infinity_method='series'):
        if not 0 <= abs(a) < 1 or mu < 0 or ell < abs(m):
            raise ValueError('Invalid Kerr or angular parameters')
        # An optional real analytic mass-squared continuation supports the
        # retarded resolvent derivative used to construct trace-driven kappa.
        # Ordinary massive-field evolution leaves this argument unset.
        msq=mu*mu if mass_squared is None else float(mass_squared)
        mu=np.sqrt(msq+0j)
        if abs(omega*omega-mu*mu) < 1e-10:
            raise ValueError('k=0 requires separate threshold asymptotics')
        self.a, self.mu, self.omega, self.m = a, mu, omega, m
        self.lam = angular_eigenvalue(ell, m, a*a*(omega*omega-mu*mu))
        self.rp = 1+np.sqrt(1-a*a)
        self.rmin, self.rmax = self.rp+offset, rmax
        self.oh = a/(2*self.rp)
        if offset <= 0 or rmax <= self.rmin:
            raise ValueError('Invalid radial interval')
        def rhs(r, state):
            d, dp, v = radial_coefficients(r, a, mu, omega, m, self.lam)
            return [state[1], -(dp*state[1]+v*state[0])/d]
        self.rhs = rhs
        r = self.rmin
        rin,din = horizon_series(r,a,mu,omega,m,self.lam)
        self.insol = solve_ivp(rhs, (r, rmax), [rin, din], dense_output=True,
                               method='DOP853', rtol=rtol, atol=rtol*1e-3)
        self.propagating = omega*omega > msq
        self.k = (np.sign(omega)*np.sqrt(omega*omega-msq) if self.propagating
                  else 1j*np.sqrt(msq-omega*omega))
        r = rmax
        k = self.k
        # A bound Up solution may be freely normalized; no infinity flux exists.
        p,f,df = infinity_series(r,a,mu,omega,m,self.lam,k)
        _,f5,_ = infinity_series(r,a,mu,omega,m,self.lam,k,order=5)
        self.series_last_term_relative=float(abs(f-f5)/max(abs(f),1e-300))
        self.infinity_method=infinity_method
        if infinity_method=='series':
            up = (np.exp(1j*k*r+p*np.log(r)-2j*k*np.log(2))*f
                  if self.propagating else 1.+0j)
            dup = up*(1j*k+p/r+df/f)
        elif infinity_method=='coulomb':
            up,dup=coulomb_boundary(r,a,mu,omega,self.lam,k,self.propagating)
        else:
            raise ValueError('Unknown infinity boundary method')
        self.upsol = solve_ivp(rhs, (rmax, self.rmin), [up, dup], dense_output=True,
                               method='DOP853', rtol=rtol, atol=rtol*1e-3)
        if not self.insol.success or not self.upsol.success:
            raise RuntimeError('Homogeneous radial integration failed')
        match = np.sqrt(self.rmin*rmax)
        self.w0 = self.wronskian(np.asarray([match]))[0]

    def wronskian(self, r):
        u, du = self.insol.sol(r)
        v, dv = self.upsol.sol(r)
        return (r*r-2*r+self.a*self.a)*(u*dv-v*du)

    def solve(self, r, source):
        """Solve for J sampled over the full finite integration interval.

        The caller must converge source-tail truncation and quadrature. Z_h
        corresponds to unit e^{-i(omega-m OmegaH)r*}; propagating Z_inf to
        e^{ik r*} r^{i mu²/k}/sqrt(r²+a²), up to boundary truncation errors.
        """
        r, source = np.asarray(r), np.asarray(source, complex)
        if (r.ndim != 1 or source.shape != r.shape or np.any(np.diff(r) <= 0)
                or r[0] != self.rmin or r[-1] != self.rmax):
            raise ValueError('Need increasing full-domain grid and equally shaped source')
        u, du = self.insol.sol(r)
        v, dv = self.upsol.sol(r)
        cu = complex_cumulative(u*source/self.w0, r)
        ci = complex_cumulative(v*source/self.w0, r)
        ci = ci[-1]-ci
        field, derivative = v*cu+u*ci, dv*cu+du*ci
        return dict(r=r, field=field, derivative=derivative, z_h=ci[0],
                    z_inf=cu[-1] if self.propagating else 0j,
                    up_coefficient=cu[-1],
                    wronskian_relative_spread=float(np.max(np.abs(self.wronskian(r)/self.w0-1))))
