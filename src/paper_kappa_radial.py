"""Direct 2023 diagonal kappa boundary-value problem, without mass derivatives.

Integrate (Delta kappa')' + V kappa = (r^2+a^2 gamma) h/2 on
both sides and impose the prescribed kappa jumps. Endpoint particular series
have their free homogeneous coefficient set to zero, then restored by matching.
This diagnostic is independent of the compact-current / mass-variation solver.
"""
import numpy as np
from numpy.polynomial import Polynomial as P
from scipy.integrate import solve_ivp
from environment_angular_diagnostic import DenseRealHarmonic
from paper_sourced_matching import orbit


def polynomials(a,w,m,lam,center=0.):
    r=P([center,1.]);d=r*r-2*r+a*a
    return d*d,d*d.deriv(),((r*r+a*a)*w-a*m)**2-d*lam,d


def infinity_data(r,a,w,m,lam,c,order=10):
    d2,dd,v,d=polynomials(a,w,m,lam)
    def coefficient(power,n,target):
        out=0j
        for degree,z in enumerate(d2.coef):
            for shift,f in ((0,-w*w),(1,2j*w*(power-n)),(2,(power-n)*(power-n-1))):
                if degree-n-shift==target:out+=z*f
        for degree,z in enumerate(dd.coef):
            if degree-n==target:out+=z*1j*w
            if degree-n-1==target:out+=z*(power-n)
        for degree,z in enumerate(v.coef):
            if degree-n==target:out+=z
        return out
    ph=-1+2j*w;pk=2j*w;h=[1.+0j]
    for n in range(1,order+2):
        target=3-n
        h.append(-sum(coefficient(ph,j,target)*z for j,z in enumerate(h))/coefficient(ph,n,target))
    source=d*P([c/2,0.,.5]);k=[]
    for n in range(order+1):
        target=3-n
        rhs=sum(z*h[degree-1-target] for degree,z in enumerate(source.coef) if 0<=degree-1-target<len(h))
        rhs-=sum(coefficient(pk,j,target)*z for j,z in enumerate(k))
        if n==1:
            if abs(rhs)>1e-7*max(1,abs(k[0])):raise RuntimeError('Infinity resonance compatibility failed')
            k.append(0j)  # Free homogeneous outgoing coefficient.
        else:k.append(rhs/coefficient(pk,n,target))
    def evaluate(coeff,power):
        pref=np.exp(1j*w*r+power*np.log(r)-2j*w*np.log(2))
        f=sum(z/r**n for n,z in enumerate(coeff));df=sum(-n*z/r**(n+1) for n,z in enumerate(coeff))
        return [pref*f,pref*((1j*w+power/r)*f+df)]
    return np.r_[evaluate(h,ph),evaluate(k,pk)]


def horizon_data(r,a,w,m,lam,c,order=8):
    rp=1+np.sqrt(1-a*a);rm=2-rp;distance=r-rp
    d2,dd,v,d=polynomials(a,w,m,lam,rp)
    # Exact horizon polynomial constant zero avoids a roundoff indicial defect.
    d=P([0.,rp-rm,1.]);d2=d*d;dd=d*d.deriv()
    rr=P([rp,1.]);v=((rr*rr+a*a)*w-a*m)**2-d*lam
    s=-1j*(2*rp*w-a*m)/(rp-rm)
    def coefficient(n,target):
        total=0j
        for degree,z in enumerate(d2.coef):
            if n-2+degree==target:total+=z*(s+n)*(s+n-1)
        for degree,z in enumerate(dd.coef):
            if n-1+degree==target:total+=z*(s+n)
        for degree,z in enumerate(v.coef):
            if n+degree==target:total+=z
        return total
    h=[1.+0j];k=[0j];source=d*(rr*rr+c)/2
    for n in range(1,order+1):
        diagonal=coefficient(n,n)
        h.append(-sum(z*coefficient(j,n) for j,z in enumerate(h))/diagonal)
        rhs=sum(z*h[n-degree] for degree,z in enumerate(source.coef) if 0<=n-degree<len(h))
        k.append((rhs-sum(z*coefficient(j,n) for j,z in enumerate(k)))/diagonal)
    def evaluate(coeff):
        pref=np.exp(s*np.log(distance));f=sum(z*distance**n for n,z in enumerate(coeff))
        df=sum(n*z*distance**(n-1) for n,z in enumerate(coeff) if n)
        return [pref*f,pref*(s/distance*f+df)]
    return np.r_[evaluate(h),evaluate(k)]


class PaperKappaRadial:
    def __init__(self,r0,a,m,ell,jumps,*,rmax=4000.,offset=1e-5,order=10,rtol=2e-12):
        if not (0<=a<1 and m>=1 and ell>=m and int(m)==m and int(ell)==ell and order>=4 and offset>0):
            raise ValueError('Diagnostic requires subextremal Kerr, positive m, ell>=m and order>=4')
        if not 1+np.sqrt(1-a*a)+offset<r0<rmax:
            raise ValueError('Orbit must lie inside both radial boundaries')
        self.r0,self.a,self.m,self.ell=r0,a,m,ell;self.w=m/(r0**1.5+a)
        harmonic=DenseRealHarmonic(0,ell,m,a*self.w)
        self.lam=harmonic.eigenvalue
        x,weights=np.polynomial.legendre.leggauss(48)
        c=a*a*2*np.pi*np.dot(weights*x*x,harmonic(np.arccos(x))**2)
        self.c=c;self.rp=1+np.sqrt(1-a*a);self.rm=2-self.rp
        self.rmin=self.rp+offset;self.rmax=rmax
        def rhs(r,y):
            d=(r-self.rp)*(r-self.rm);K=(r*r+a*a)*self.w-a*m;V=K*K/d-self.lam
            h,hp,k,kp=y
            return np.array([hp,(-2*(r-1)*hp-V*h)/d,kp,((r*r+c)*h/2-2*(r-1)*kp-V*k)/d])
        self.rhs=rhs
        args=dict(method='DOP853',rtol=rtol,atol=rtol*1e-4,dense_output=True)
        self.inner=solve_ivp(rhs,(self.rmin,r0),horizon_data(self.rmin,a,self.w,m,self.lam,c,order),**args)
        self.outer=solve_ivp(rhs,(rmax,r0),infinity_data(rmax,a,self.w,m,self.lam,c,order),**args)
        if not self.inner.success or not self.outer.success:raise RuntimeError('Direct kappa integration failed')
        inside=self.inner.y[:,-1];outside=self.outer.y[:,-1]
        A=np.column_stack([outside[:2],-inside[:2]])
        _,ut,_=orbit(r0,a);delta=(r0-self.rp)*(r0-self.rm)
        self.hamps=np.linalg.solve(A,[0.,-16*np.pi*harmonic(np.pi/2)/(ut*delta)])
        known=outside[2:]*self.hamps[0]-inside[2:]*self.hamps[1]
        self.kamps=np.linalg.solve(A,np.asarray(jumps)-known)
        self.jump_residual=A@self.kamps+known-jumps

    def state(self,r):
        idx=1 if r<self.r0 else 0
        sol=self.inner if idx else self.outer;y=sol.sol(r)
        h,hp=y[:2]*self.hamps[idx]
        k,kp=y[2:]*self.hamps[idx]+y[:2]*self.kamps[idx]
        return np.array([h,hp,k,kp])
