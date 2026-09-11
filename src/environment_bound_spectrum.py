"""Complex-frequency shooting of low-lying massive scalar Kerr bound modes.

Use continuation to identify a branch. A small log-derivative residual alone
does not establish the overtone label or the accuracy of a tiny growth rate.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from environment_radial import horizon_series,infinity_series


def complex_angular_eigenvalue(ell,m,c2,size=24):
    ls=np.arange(abs(m),abs(m)+size+1)
    C=np.zeros((size+1,size+1))
    for j,l in enumerate(ls[:-1]):
        C[j,j+1]=C[j+1,j]=np.sqrt(((l+1)**2-m*m)/((2*l+1)*(2*l+3)))
    matrix=np.diag(ls[:size]*(ls[:size]+1))-c2*(C@C)[:size,:size]
    return min(np.linalg.eigvals(matrix),key=lambda z:abs(z-ell*(ell+1)))


def quasibound_mode(a,mu,ell,m,*,seed=None,outer_efolds=30.,offset=1e-5,rtol=2e-10):
    if not 0<=a<1 or not .09<=mu<=.31 or not abs(m)<=ell<=2:
        raise ValueError('Targeted range: alpha=.09..31, ell<=2, |m|<=ell, 0<=a<1')
    rp=1+np.sqrt(1-a*a);rm=1-np.sqrt(1-a*a);rmin=rp+offset
    match=(ell+1)/mu**2;rmax=outer_efolds*(ell+1)/mu**2
    if not rp<rmin<match<rmax:raise ValueError('Invalid shooting interval')
    def residual(omega):
        lam=complex_angular_eigenvalue(ell,m,a*a*(omega*omega-mu*mu))
        k=1j*np.sqrt(mu*mu-omega*omega+0j)
        if k.imag<=0:raise ValueError('Infinity branch must decay spatially')
        def rhs(r,y):
            # Factored Delta avoids cancellation of r^2-2r+a^2 close to H.
            d=(r-rp)*(r-rm);dp=2*(r-1)
            K=omega*(r-rp)*(r+rp)+(2*rp*omega-a*m)
            v=K*K/d-mu*mu*r*r-a*a*omega*omega+2*a*m*omega-lam
            return [y[1],-(dp*y[1]+v*y[0])/d]
        options=dict(method='DOP853',rtol=rtol,atol=rtol*1e-3)
        def inner_rhs(x,y):
            distance=np.exp(x);r=rp+distance;b=rp-rm
            K=omega*distance*(r+rp)+(2*rp*omega-a*m)
            A=mu*mu*r*r+a*a*omega*omega-2*a*m*omega+lam
            return [y[1],-distance/(distance+b)*y[1]-(K*K/(distance+b)**2-distance*A/(distance+b))*y[0]]
        rin,din=horizon_series(rmin,a,mu,omega,m,lam)
        left=solve_ivp(inner_rhs,(np.log(rmin-rp),np.log(match-rp)),[rin,(rmin-rp)*din],**options)
        p,f,df=infinity_series(rmax,a,mu,omega,m,lam,k)
        right=solve_ivp(rhs,(rmax,match),[1.+0j,1j*k+p/rmax+df/f],**options)
        if not left.success or not right.success:raise RuntimeError('Radial shooting failed')
        return left.y[1,-1]/((match-rp)*left.y[0,-1])-right.y[1,-1]/right.y[0,-1]
    if seed is None:
        wr=mu*(1-mu*mu/(2*(ell+1)**2))
        seed=wr-1j*np.sign(wr-m*a/(2*rp))*mu**(4*ell+6)
    def equations(x):
        value=residual(complex(*x));return [value.real,value.imag]
    answer=root(equations,[seed.real,seed.imag],tol=1e-11)
    omega=complex(*answer.x)
    for polish in range(6):
        value=residual(omega)
        if abs(value)<1e-13:break
        step=1e-6
        derivative=(residual(omega+step)-residual(omega-step))/(2*step)
        omega-=value/derivative
    value=residual(omega)
    if not 0<omega.real<mu or abs(value)>1e-10:
        raise RuntimeError(f'Unresolved root {omega}, residual {value}')
    return dict(a=a,alpha=mu,ell=ell,m=m,omega=[omega.real,omega.imag],
        residual=[value.real,value.imag],match_radius=match,rmax=rmax,
        horizon_offset=offset,rtol=rtol,outer_efolds=outer_efolds,
        root_solver_success=bool(answer.success),newton_polish_iterations=polish,
        branch_label='Not determined by a single solve; use parameter continuation')


def continue_fundamental(a,mu,ell,m):
    masses=list(np.arange(.1,mu-1e-9,.05))+[mu]
    rows=[];seed=None;previous=None
    for mass in masses:
        if rows:
            old=complex(*rows[-1]['omega'])
            seed=mass+(old-previous)*(mass/previous)**3
        row=quasibound_mode(a,float(mass),ell,m,seed=seed)
        rows.append(row);previous=mass
    return rows
