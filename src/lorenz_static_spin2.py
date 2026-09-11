"""Static vacuum Hertz-to-Lorenz construction from 2306.16459v3 Sec. IV.

Legendre P/Q Hertz branches, with optional point-particle Weyl normalization.
Chi has zero data at a reference radius. Homogeneous matching and physical
boundary adaptation are not supplied. The complex vacuum metric is made real
by adding its conjugate with the calibrated Hertz amplitude.
"""
from functools import lru_cache
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import eval_legendre,lqmn
from numpy.polynomial import Legendre,Polynomial
from lorenz_jet import Jet
from lorenz_ghp import KerrGHP
from lorenz_tensor import tensor_divergence,vector_covariant_derivative


def _legendre_jet(x,ell,derivative=0):
    coefficients=Legendre.basis(ell).convert(kind=Polynomial).deriv(derivative).coef
    result=Jet(0.,x.order)
    for c in coefficients[::-1]:
        result=result*x+c
    return result


def _angular(g,ell,spin=-2):
    norm=np.sqrt((2*ell+1)/(4*np.pi))
    if spin==0:return norm*_legendre_jet(g.theta.cos(),ell)
    norm*=np.sqrt(math.factorial(ell-2)/math.factorial(ell+2))
    return norm*g.theta.sin()**2*_legendre_jet(g.theta.cos(),ell,2)


def _radial(r,a,ell,branch):
    b=np.sqrt(1-a*a)
    x=(r-1)/b
    lam=ell*(ell+1)
    if branch=='P':
        R=eval_legendre(ell,x)
        Rp=ell*(x*R-eval_legendre(ell-1,x))/(x*x-1)/b
    elif branch=='Q':
        values,derivatives=lqmn(0,ell,x)
        R,Rp=values[0,ell],derivatives[0,ell]/b
    else:raise ValueError('Static radial branch must be P or Q')
    delta=r*r-2*r+a*a
    return delta*(lam*R-2*(r-1)*Rp),(lam-2)*delta*Rp,lam*(lam-2)*R


def _radial_jet(g,ell,branch):
    p,dp,_=_radial(float(g.r.value.real),g.a,ell,branch)
    result=Jet(p,g.order); result.c[1,0]=dp
    for n in range(g.order-1):
        rhs=(2*(g.r-1)*result.derivative(0)+(ell*(ell+1)-2)*result)/g.delta
        result.c[n+2,0]=rhs.c[n,0]/((n+1)*(n+2))
    return result


@lru_cache(maxsize=32)
def _couplings(ell):
    x,w=np.polynomial.legendre.leggauss(max(20,ell+5))
    values=[]
    for t in np.arccos(x):
        g=KerrGHP(6.,t,0.,order=3)
        s,c=g.theta.sin(),g.theta.cos()
        Y=_angular(g,ell)
        L2=Y.derivative(1)+2*c/s*Y
        L12=L2.derivative(1)+c/s*L2
        values.append([L12.value,(s*c*L2+2*s*s*Y).value,(c*L12).value,(s*L2).value])
    values=np.asarray(values)
    return tuple((j,tuple(2*np.pi*np.einsum('i,ij->j',w*np.sqrt((2*j+1)/(4*np.pi))*eval_legendre(j,x),values)))
                 for j in range(max(0,ell-2),ell+3))


@lru_cache(maxsize=128)
def _chi_state(r,a,ell,branch,reference_radius):
    modes=_couplings(ell); lam=ell*(ell+1)
    if r==reference_radius:return np.zeros(2*len(modes),complex)
    def rhs(x,y):
        p,dp,ddp=_radial(x,a,ell,branch)
        delta=x*x-2*x+a*a
        values=[]
        for n,(j,(A,B,C,D)) in enumerate(modes):
            source=(-A*(x*dp-2*p)+a*a*B*ddp+1j*a*(C*dp+D*x*ddp))/lam
            values.extend((y[2*n+1],(j*(j+1)*y[2*n]+source-2*(x-1)*y[2*n+1])/delta))
        return values
    solution=solve_ivp(rhs,(reference_radius,r),np.zeros(2*len(modes),complex),method='DOP853',rtol=2e-12,atol=2e-14)
    if not solution.success:raise RuntimeError(solution.message)
    return solution.y[:,-1]


def static_spin2_metric(r,theta,a=.6,ell=2,branch='P',order=8,reference_radius=6.,radiation_gauge=False):
    if ell<2 or not 0<=a<1 or min(r,reference_radius)<=1+np.sqrt(1-a*a):
        raise ValueError('Static ell>=2 exterior vacuum data required')
    g=KerrGHP(r,theta,a,omega=0.,m=0,order=order)
    p=_radial_jet(g,ell,branch)
    Y=_angular(g,ell)
    psi=p*Y
    s,c=g.theta.sin(),g.theta.cos()
    rho=g.r+1j*a*c; rhoc=g.r-1j*a*c
    L2=psi.derivative(1)+2*c/s*psi
    A=L2.derivative(1)+c/s*L2-2j*a*s/rhoc*L2
    B=psi.derivative(0).derivative(0)-2/rhoc*psi.derivative(0)
    C=L2.derivative(0)-(L2-1j*a*s*psi.derivative(0))/rho-(L2+1j*a*s*psi.derivative(0))/rhoc
    l=g.tetrad[0]
    m=[1j*a*s,Jet(0.,order),Jet(1.,order),1j/s]
    upper=[[-(l[i]*l[j]*A+m[i]*m[j]*B-(l[i]*m[j]+m[i]*l[j])*C)/(2*rho*rho)
            for j in range(4)] for i in range(4)]
    h=[[sum(g.g[i][u]*g.g[j][v]*upper[u][v] for u in range(4) for v in range(4))
        for j in range(4)] for i in range(4)]
    if radiation_gauge:return g,h
    lam=ell*(ell+1)
    # Minus sign is fixed by the sourced spin-1 equation for rhoc*HU,
    # not by fitting a metric residual; see STATIC_SPIN2.md for an explicit
    # Schwarzschild polynomial substitution and the printed-sign discrepancy.
    HU=-(rhoc*L2.derivative(0)-2*(L2+1j*a*s*psi.derivative(0)))/(2*lam*rhoc)
    twoform=[[sum(g.g[i][u]*g.g[j][v]*(l[u]*m[v]-m[u]*l[v])*HU/g.sigma
                  for u in range(4) for v in range(4)) for j in range(4)] for i in range(4)]
    chi=Jet(0.,order)
    state=_chi_state(r,a,ell,branch,reference_radius)
    for n,(j,(aa,bb,cc,dd)) in enumerate(_couplings(ell)):
        radial=Jet(state[2*n],order); radial.c[1,0]=state[2*n+1]
        source=(-aa*(g.r*p.derivative(0)-2*p)+a*a*bb*p.derivative(0).derivative(0)
                +1j*a*(cc*p.derivative(0)+dd*g.r*p.derivative(0).derivative(0)))/lam
        for k in range(order-1):
            rhs=(j*(j+1)*radial+source-2*(g.r-1)*radial.derivative(0))/g.delta
            radial.c[k+2,0]=rhs.c[k,0]/((k+1)*(k+2))
        chi+=radial*_angular(g,j,spin=0)
    div=tensor_divergence(g,twoform)
    xi=[-rhoc*rhoc*div[i]-g.partial(chi,i) for i in range(4)]
    derivative=vector_covariant_derivative(g,xi)
    return g,[[h[i][j]-derivative[i][j]-derivative[j][i] for j in range(4)] for i in range(4)]


@lru_cache(maxsize=128)
def static_hertz_amplitude(r0,a,ell,branch):
    """Invert psi0=(1/2) D^4 conjugate(psi_Hertz), m=omega=0.

    The P/Q radial normalization is arbitrary, so compare it to the pinned
    Teukolsky normalization at one vacuum radius. This fixes curvature only;
    it does not fix completion or the homogeneous scalar gauge freedom.
    """
    from pybhpt.radial import RadialTeukolsky
    from lorenz_weyl import weyl_amplitudes
    rp=1+np.sqrt(1-a*a)
    sample=(rp+r0)/2 if branch=='P' else 2*r0
    bc,index=('In',1) if branch=='P' else ('Up',0)
    radial=RadialTeukolsky(2,ell,0,a,0.,np.array([sample]))
    radial.solve(bc=bc)
    target=weyl_amplitudes(r0,a,ell,0)[2][index]*radial.radialsolution(bc,0)
    p=_radial(sample,a,ell,branch)[0]
    delta=sample*sample-2*sample+a*a
    lam=ell*(ell+1)
    return (2*target/(lam*(lam-2)*p/delta**2)).conjugate()


def sourced_static_spin2(r,theta,r0=6.,a=.6,ell=2,order=8,circular_symmetry=True):
    """Curvature-normalized circular static piece, not a matched metric.

    Average with the Kerr isometry (t,phi)->(-t,-phi). The stationary
    circular source is invariant under this pullback. A single IRG seed
    otherwise retains a curvature-free, symmetry-odd gauge contribution.
    """
    if r==r0:raise ValueError('Use a one-sided vacuum point')
    branch='P' if r<r0 else 'Q'
    amplitude=static_hertz_amplitude(r0,a,ell,branch)
    g,h=static_spin2_metric(r,theta,a=a,ell=ell,branch=branch,order=order,reference_radius=r0)
    real=[[amplitude*v+amplitude.conjugate()*v.conjugate() for v in row] for row in h]
    if circular_symmetry:
        parity=(-1,1,1,-1)
        real=[[v*(1+parity[i]*parity[j])/2 for j,v in enumerate(row)] for i,row in enumerate(real)]
    return g,real
