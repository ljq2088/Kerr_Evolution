"""In/Up Green integrals for the static scalar and r-theta two-form ODEs.

These fix the auxiliary radial boundary freedom, not the metric completion.
Integration errors reported here exclude source and floating-point errors.
"""
from functools import lru_cache
import numpy as np
from scipy.integrate import quad_vec
from scipy.special import eval_legendre,lqmn


class StaticGreen:
    def __init__(self,a,ell,source,pivot=6.,kind='scalar',breaks=(),rtol=2e-11,atol=2e-12):
        if not 0<=a<1 or int(ell)!=ell or ell<0 or kind not in ('scalar','twoform'):
            raise ValueError('Subextremal Kerr, integer ell>=0, scalar/twoform required')
        if kind=='twoform' and ell==0:
            raise ValueError('There is no two-form ell=0 mode')
        self.a,self.ell,self.source,self.kind=a,ell,source,kind
        self.b=np.sqrt(1-a*a);self.rp=1+self.b
        if pivot<=self.rp:
            raise ValueError('Normalization pivot must be outside the horizon')
        self.pivot=pivot;self.breaks=tuple(sorted(set(breaks)))
        self.rtol,self.atol=rtol,atol
        self.polynomial=np.polynomial.legendre.Legendre.basis(ell)
        raw=self._raw(pivot)
        self.pnorm,self.qnorm=raw[0],raw[2]
        self.wronskian=(-self.b if kind=='scalar' else ell*(ell+1)*self.b)/(self.pnorm*self.qnorm)

    def _raw(self,r):
        x=(r-1)/self.b
        P=eval_legendre(self.ell,x)
        dP=self.polynomial.deriv()(x)/self.b
        q,dq=lqmn(0,self.ell,x)
        Q,dQ=q[0,self.ell],dq[0,self.ell]/self.b
        if self.kind=='twoform':
            delta=r*r-2*r+self.a*self.a
            lam=self.ell*(self.ell+1)
            return np.array([delta*dP,lam*P,delta*dQ,lam*Q])
        return np.array([P,dP,Q,dQ])

    def basis(self,r):
        return self._raw(r)/np.array([self.pnorm,self.pnorm,self.qnorm,self.qnorm])

    def _integral(self,lo,hi,branch):
        if lo==hi:return 0j,0.
        points=[lo]+[v for v in self.breaks if lo<v<hi]+[hi]
        value,error=0j,0.
        def integrand(r):
            source=self.source(r)
            if self.kind=='twoform':source/=r*r-2*r+self.a*self.a
            result=self.basis(r)[branch]*source
            if not np.isfinite(result):
                raise FloatingPointError(f'Nonfinite static Green integrand at r={r}')
            return result
        for left,right in zip(points[:-1],points[1:]):
            v,e,info=quad_vec(integrand,left,right,epsrel=self.rtol,epsabs=self.atol,
                              limit=500,full_output=True)
            if not info.success:
                raise RuntimeError(f'Static Green quadrature failed: {info.message}')
            value+=v;error+=e
        return value,error

    @lru_cache(maxsize=512)
    def evaluate(self,r):
        if r<=self.rp:
            raise ValueError('Evaluate at an exterior vacuum point')
        incoming,e_in=self._integral(self.rp,r,0)
        outgoing,e_up=self._integral(r,np.inf,2)
        P,dP,Q,dQ=self.basis(r)
        R=(Q*incoming+P*outgoing)/self.wronskian
        dR=(dQ*incoming+dP*outgoing)/self.wronskian
        delta=r*r-2*r+self.a*self.a
        ddR=(self.ell*(self.ell+1)*R+self.source(r)
             -(2*(r-1)*dR if self.kind=='scalar' else 0))/delta
        error=(abs(Q)*e_in+abs(P)*e_up)/abs(self.wronskian)
        return np.array([R,dR,ddR]),float(error)


@lru_cache(maxsize=128)
def _trace_green(r0,a,ell,j,kind):
    from lorenz_static_trace import static_trace_radial
    from lorenz_static_gauge import angular_couplings
    _,c,d=next(v for v in angular_couplings(ell) if v[0]==j)
    def source(r):
        radial,_=static_trace_radial(r,r0,a,ell,side='outside' if r>=r0 else 'inside')
        if kind=='scalar':
            return radial[0]*((r*r if j==ell else 0)+a*a*c)/2
        return (r*r-2*r+a*a)*radial[1]*((r*r if j==ell else 0)+a*a*d)/(2*ell*(ell+1))
    return StaticGreen(a,j,source,pivot=r0,kind=kind,breaks=(r0,))


def regular_trace_state(r,r0,a,ell):
    from lorenz_static_gauge import angular_couplings
    state=[]
    for j,_,_ in angular_couplings(ell):
        k,_=_trace_green(r0,a,ell,j,'scalar').evaluate(r)
        b=_trace_green(r0,a,ell,j,'twoform').evaluate(r)[0] if j else np.zeros(3)
        state.extend((k[0],k[1],b[0],b[1]))
    return np.asarray(state)


@lru_cache(maxsize=128)
def _chi_green(r0,a,ell,j):
    from lorenz_static_spin2 import _couplings,_radial,static_hertz_amplitude
    _,(aa,bb,cc,dd)=next(v for v in _couplings(ell) if v[0]==j)
    amplitudes={branch:static_hertz_amplitude(r0,a,ell,branch) for branch in ('P','Q')}
    def source(r):
        branch='P' if r<r0 else 'Q'
        p,dp,ddp=np.asarray(_radial(r,a,ell,branch))*amplitudes[branch]
        return (-aa*(r*dp-2*p)+a*a*bb*ddp+1j*a*(cc*dp+dd*r*ddp))/(ell*(ell+1))
    return StaticGreen(a,j,source,pivot=r0,breaks=(r0,))


def regular_chi_state(r,r0,a,ell):
    from lorenz_static_spin2 import _couplings
    return np.asarray([value for j,_ in _couplings(ell) for value in _chi_green(r0,a,ell,j).evaluate(r)[0][:2]])


@lru_cache(maxsize=128)
def _scalar_basis(a,ell,r0):
    return StaticGreen(a,ell,lambda r:0.,pivot=r0)


def free_scalar_amplitudes(a,ell,r0,jump_value,jump_derivative):
    """Horizon-regular P / infinity-decaying Q amplitudes for prescribed jumps.

    Both radial basis functions equal one at r0. Jumps are Up minus In.
    """
    _,dp,_,dq=_scalar_basis(a,ell,r0).basis(r0)
    incoming=(jump_derivative-dq*jump_value)/(dq-dp)
    return incoming,incoming+jump_value


def matched_free_scalar_metric(r,theta,r0,a,ell,jump_value,jump_derivative,order=6,side=None):
    from lorenz_ghp import KerrGHP
    from lorenz_jet import Jet
    from lorenz_static_gauge import _angular_jet
    from lorenz_tensor import vector_covariant_derivative
    if ell<2 or ell%2 or (r==r0 and side not in ('inside','outside')):
        raise ValueError('Even ell>=2 required; select a side if r==r0')
    amplitudes=free_scalar_amplitudes(a,ell,r0,jump_value,jump_derivative)
    inside=r<r0 or (r==r0 and side=='inside')
    basis=_scalar_basis(a,ell,r0).basis(r)
    index=0 if inside else 2
    R,dR=basis[index:index+2]*amplitudes[0 if inside else 1]
    g=KerrGHP(r,theta,a,order=order)
    radial=Jet(R,order);radial.c[1,0]=dR
    for n in range(order-1):
        rhs=(ell*(ell+1)*radial-2*(g.r-1)*radial.derivative(0))/g.delta
        radial.c[n+2,0]=rhs.c[n,0]/((n+1)*(n+2))
    field=radial*_angular_jet(g,ell)
    h=vector_covariant_derivative(g,[g.partial(field,j) for j in range(4)])
    return g,h,field


@lru_cache(maxsize=16)
def _completion_y2_green(a):
    b=np.sqrt(1-a*a);rp,rm=1+b,1-b
    def source(r):
        return 2*a*a/3*np.log((r-rp)/(r-rm))/b
    return StaticGreen(a,2,source,pivot=6.)


def regular_completion_y2(a,r):
    """Coefficient of Legendre P2 in y; excludes the growing homogeneous r² term."""
    return _completion_y2_green(a).evaluate(r)[0]
