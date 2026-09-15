"""Independent normalized ingoing massive-scalar radial series on Kerr.

No production radial coefficients, boundary expansions or ODE integration
are called. Uses the factored KG equation and Frobenius recurrence in
x=(r-r+)/(r-r-). The phase normalization is fixed analytically, never fitted.
"""
import mpmath as mp
import numpy as np


class LeaverIngoing:
    def __init__(self,a,mu,omega,ell,m,terms=400,dps=70):
        self.dps=dps
        with mp.workdps(dps):
            a,mu,w=[mp.mpf(str(v)) for v in (a,mu,omega)]
            rp=1+mp.sqrt(1-a*a);rm=2-rp;d=rp-rm
            q=-mp.sqrt(mu*mu-w*w) if abs(w)<mu else 1j*mp.sign(w)*mp.sqrt(w*w-mu*mu)
            beta=(mu*mu-2*w*w)/q-1;s=-1j*(2*rp*w-a*m)/d
            lmin=abs(m)+(ell-abs(m))%2;size=max(12,(ell-lmin)//2+8)
            C=lambda l:mp.sqrt(mp.mpf(l*l-m*m)/(4*l*l-1)) if l>abs(m) else mp.mpf(0)
            matrix=mp.matrix(size);c2=a*a*(w*w-mu*mu)
            for i in range(size):
                l=lmin+2*i;matrix[i,i]=l*(l+1)-c2*(C(l+1)**2+C(l)**2)
                if i<size-1:matrix[i,i+1]=matrix[i+1,i]=-c2*C(l+1)*C(l+2)
            angular=mp.eigsy(matrix,eigvals_only=True)[(ell-lmin)//2]
            sep=angular+a*a*w*w-2*a*m*w
            def c(x):
                r=(rp-rm*x)/(1-x);delta=(r-rp)*(r-rm);K=(r*r+a*a)*w-a*m
                lp=s/x+beta/(1-x)+q*d/(1-x)**2
                lpp=-s/x**2+beta/(1-x)**2+2*q*d/(1-x)**3
                return (1-x)**2*(x*(lpp+lp*lp)+lp)+K*K/delta-mu*mu*r*r-sep
            x1,x2=mp.mpf('.25'),mp.mpf('.75');c1=(c(x2)-c(x1))/(x2-x1);c0=c(x1)-c1*x1
            if abs(c(mp.mpf('.5'))-c0-c1/2)>mp.mpf(10)**(-dps+10):raise RuntimeError('Radial reduction is not linear')
            b0=1+2*s;b1=-2*b0+2*(q*d+beta);b2=b0-2*beta
            coeff=[mp.mpc(1)]
            for n in range(terms):
                al=(n+1)*(n+b0);be=-2*n*(n-1)+b1*n+c0;ga=(n-1)*(n-2)+b2*(n-1)+c1
                coeff.append(-(be*coeff[-1]+(ga*coeff[-2] if n else 0))/al)
            dw=w-a*m/(2*rp);constant=rp-2*rp/d*mp.log(2)-2*rm/d*mp.log(d/2)
            normalization=mp.exp(-q*rp-1j*dw*constant)*d**(s-beta)
            self.rp,self.rm,self.q,self.beta,self.s,self.norm=rp,rm,q,beta,s,normalization
            self.coeff=coeff;self.angular=angular

    def state_mp(self,r):
        with mp.workdps(self.dps):
            r=mp.mpf(str(r));x=(r-self.rp)/(r-self.rm)
            if not 0<x<1:raise ValueError('Series requires an exterior radius')
            f=mp.polyval(list(reversed(self.coeff)),x)
            fx=mp.polyval([n*self.coeff[n] for n in range(len(self.coeff)-1,0,-1)],x)
            xp=(self.rp-self.rm)/(r-self.rm)**2
            pref=self.norm*mp.exp(self.q*r)*(r-self.rm)**self.beta*x**self.s
            return pref*f,pref*(fx*xp+f*(self.q+self.beta/(r-self.rm)+self.s*xp/x))

    def state(self,r):
        return np.array([complex(v) for v in self.state_mp(r)])
