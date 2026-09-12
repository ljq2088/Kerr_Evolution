"""Differentiated real-frequency scalar radial solutions, diagnostic only.

Auxiliary mass-squared derivatives are propagated with the original ODE.
Both finite series boundary data are differentiated by forward arithmetic.
Production Lorenz reconstruction continues to use its existing implementation.
"""
import cmath
import math
import numpy as np
from scipy.integrate import solve_ivp
from environment_angular_variation import angular_eigenpair_mass_derivative


class _Dual:
    def __init__(self,value,derivative=0):
        self.value=complex(value);self.derivative=complex(derivative)
    @staticmethod
    def cast(other):return other if isinstance(other,_Dual) else _Dual(other)
    def __add__(self,other):
        o=self.cast(other);return _Dual(self.value+o.value,self.derivative+o.derivative)
    __radd__=__add__
    def __neg__(self):return _Dual(-self.value,-self.derivative)
    def __sub__(self,other):return self+-self.cast(other)
    def __rsub__(self,other):return self.cast(other)+-self
    def __mul__(self,other):
        o=self.cast(other);return _Dual(self.value*o.value,self.derivative*o.value+self.value*o.derivative)
    __rmul__=__mul__
    def __truediv__(self,other):
        o=self.cast(other);return _Dual(self.value/o.value,(self.derivative*o.value-self.value*o.derivative)/(o.value*o.value))
    def __rtruediv__(self,other):return self.cast(other)/self
    def exp(self):
        value=cmath.exp(self.value);return _Dual(value,value*self.derivative)


def _pack(value,derivative):
    return np.array([value.value,derivative.value,value.derivative,derivative.derivative],complex)


def horizon_mass_boundary(r,a,omega,m,eigenvalue,eigenvalue_mass,mass_squared=0.,order=4):
    rp=1+math.sqrt(1-a*a);rm=2-rp;d=rp-rm;x=r-rp
    if x<=0:raise ValueError('Horizon boundary must be exterior')
    nu=_Dual(mass_squared,1);lam=_Dual(eigenvalue,eigenvalue_mass)
    k0=2*rp*omega-a*m
    k2=np.polynomial.polynomial.polymul([k0,2*rp*omega,omega],[k0,2*rp*omega,omega])
    b0=nu*rp*rp+a*a*omega*omega-2*a*m*omega+lam
    b1=2*rp*nu;b2=nu
    delta_b=[_Dual(0),d*b0,b0+d*b1,b1+d*b2,b2]
    potential=[_Dual(float(v))-b for v,b in zip(k2,delta_b)]
    dw=omega-m*a/(2*rp);s=-1j*dw*2*rp/d
    def coefficient(n,target):
        value=_Dual(0)
        for j,v in ((2,d*d),(3,2*d),(4,1)):
            if n-2+j==target:value+=v*(s+n)*(s+n-1)
        for j,v in ((1,d*d),(2,3*d),(3,2)):
            if n-1+j==target:value+=v*(s+n)
        if 0<=target-n<len(potential):value+=potential[target-n]
        return value
    coeff=[_Dual(1)]
    for n in range(1,order+1):
        coeff.append(-sum(c*coefficient(j,n) for j,c in enumerate(coeff))/coefficient(n,n))
    f=sum(c*x**n for n,c in enumerate(coeff))
    df=sum(n*c*x**(n-1) for n,c in enumerate(coeff) if n)
    constant=rp-2*rp/d*math.log(2)-2*rm/d*math.log(d/2)
    phase=cmath.exp(s*math.log(x)-1j*dw*constant)
    return _pack(phase*f,phase*(s*f/x+df))


def infinity_mass_boundary(r,a,omega,m,eigenvalue,eigenvalue_mass,mass_squared=0.,order=6):
    if omega==0 or mass_squared>=omega*omega:
        raise ValueError('This variation requires the propagating auxiliary branch')
    nu=_Dual(mass_squared,1);lam=_Dual(eigenvalue,eigenvalue_mass)
    kval=math.copysign(math.sqrt(omega*omega-mass_squared),omega)
    k=_Dual(kval,-1/(2*kval));p=-1+1j*(2*omega*omega-nu)/k
    d2={4:1.,3:-4.,2:4+2*a*a,1:-4*a*a,0:a**4}
    dd={3:2.,2:-6.,1:4+2*a*a,0:-2*a*a}
    C=a*a*omega*omega-2*a*m*omega+lam;b=a*a*omega-a*m
    potential={4:k*k,3:2*nu,2:2*omega*b-C-a*a*nu,1:2*C,0:b*b-a*a*C}
    def coefficient(n,target):
        total=_Dual(0)
        for degree,v in d2.items():
            for shift,f in ((0,-k*k),(1,2j*k*(p-n)),(2,(p-n)*(p-n-1))):
                if degree-n-shift==target:total+=v*f
        for degree,v in dd.items():
            for shift,f in ((0,1j*k),(1,p-n)):
                if degree-n-shift==target:total+=v*f
        return total+potential.get(target+n,0.)
    coeff=[_Dual(1)]
    for n in range(1,order+1):
        target=3-n
        coeff.append(-sum(c*coefficient(j,target) for j,c in enumerate(coeff))/coefficient(n,target))
    f=sum(c/r**n for n,c in enumerate(coeff));df=sum(-n*c/r**(n+1) for n,c in enumerate(coeff))
    value=(1j*k*r+p*math.log(r)-2j*k*math.log(2)).exp()*f
    return _pack(value,value*(1j*k+p/r+df/f))


class RadialMassVariation:
    def __init__(self,a,omega,ell,m,*,mass_squared=0.,rmax=2000.,offset=1e-4,rtol=1e-12):
        if not 0<=abs(a)<1 or abs(omega)<1e-12 or mass_squared>=omega*omega:
            raise ValueError('Require subextremal Kerr and nonzero propagating auxiliary frequency')
        if not all(np.isfinite(v) for v in (a,omega,mass_squared,rmax,offset,rtol)) or offset<=0 or rtol<=0:
            raise ValueError('Finite positive numerical controls required')
        self.a=float(a);self.omega=float(omega);self.m=int(m);self.mass_squared=float(mass_squared)
        lam,_,dlam,_,_=angular_eigenpair_mass_derivative(ell,m,a,omega,mass_squared)
        self.lam,self.lam_mass=lam,dlam;self.rp=1+math.sqrt(1-a*a)
        self.rmin=self.rp+offset;self.rmax=rmax
        if rmax<=self.rmin:raise ValueError('Invalid radial interval')
        def rhs(r,y):
            R,P,U,Q=y;delta=r*r-2*r+a*a;dp=2*(r-1)
            K=(r*r+a*a)*omega-a*m
            V=K*K/delta-mass_squared*r*r-a*a*omega*omega+2*a*m*omega-lam
            return [P,-(dp*P+V*R)/delta,Q,(-dp*Q-V*U+(r*r+dlam)*R)/delta]
        self.rhs=rhs
        initial=horizon_mass_boundary(self.rmin,a,omega,m,lam,dlam,mass_squared)
        outer=infinity_mass_boundary(rmax,a,omega,m,lam,dlam,mass_squared)
        self.insol=solve_ivp(rhs,(self.rmin,rmax),initial,method='DOP853',rtol=rtol,atol=rtol*1e-3,dense_output=True)
        self.upsol=solve_ivp(rhs,(rmax,self.rmin),outer,method='DOP853',rtol=rtol,atol=rtol*1e-3,dense_output=True)
        if not self.insol.success or not self.upsol.success:raise RuntimeError('Variational radial integration failed')
        self.w0,self.w0_mass=self.wronskian(np.array([math.sqrt(self.rmin*rmax)]))[:,0]
    def wronskian(self,r):
        r=np.asarray(r);u,du,U,dU=self.insol.sol(r);v,dv,V,dV=self.upsol.sol(r)
        delta=r*r-2*r+self.a*self.a
        return np.array([delta*(u*dv-v*du),delta*(U*dv+u*dV-V*du-v*dU)])
