"""Independent real-frequency |211> threshold cloud from a Leaver series.

M=1, m=1. The radial equation is reduced with x=(r-r+)/(r-r-),
R=exp(-k r)*(r-r-)^beta*sum a_n*x^n. Minimal recurrence coefficients
are selected by a continued fraction, not by the production radial shooting.
"""
import mpmath as mp
import numpy as np
from scipy.special import lpmv,gammaln,roots_laguerre


class LeaverThresholdCloud:
    def __init__(self,mu=.3,terms=200,angular_terms=10,dps=50):
        if not .19<=mu<=.31 or terms<30 or angular_terms<3:raise ValueError('Validated |211> range and adequate truncations required')
        self.dps=dps;self.terms=terms;self.angular_terms=angular_terms
        with mp.workdps(dps):
            self.mu=mp.mpf(str(mu))
            spin=lambda w:4*w/(1+4*w*w)
            seeds=(spin(self.mu*(1-mp.mpf('1.2')*self.mu**2/8)),
                   spin(self.mu*(1-mp.mpf('0.8')*self.mu**2/8)))
            self.a=mp.findroot(lambda a:self.recurrence(a)[0],seeds,solver='anderson',
                               tol=mp.mpf(10)**(-dps+8),maxsteps=100)
            self.residual,data=self.recurrence(self.a)
            self.rp,self.rm,self.omega,self.k,self.beta,self.angular_coefficients,ratios=data
            coefficients=[mp.mpf(1)]
            for ratio in ratios:coefficients.append(coefficients[-1]*ratio)
            self.coefficients=coefficients

    def recurrence(self,a):
        rp=1+mp.sqrt(1-a*a);rm=2-rp;d=rp-rm;w=a/(2*rp)
        k=mp.sqrt(self.mu*self.mu-w*w);beta=(2*w*w-self.mu*self.mu)/k-1
        c2=a*a*(w*w-self.mu*self.mu)
        C=lambda l:mp.sqrt(mp.mpf(l*l-1)/(4*l*l-1)) if l>1 else mp.mpf(0)
        matrix=mp.matrix(self.angular_terms)
        for i in range(self.angular_terms):
            ell=1+2*i;matrix[i,i]=ell*(ell+1)-c2*(C(ell+1)**2+C(ell)**2)
            if i<self.angular_terms-1:
                matrix[i,i+1]=matrix[i+1,i]=-c2*C(ell+1)*C(ell+2)
        values,vectors=mp.eigsy(matrix);A=values[0]
        angular=[vectors[i,0]*mp.sign(vectors[0,0]) for i in range(self.angular_terms)]
        separation=A+a*a*w*w-2*a*w
        b1=-2+2*(beta-k*d);b2=1-2*beta
        c0=beta-k*d-self.mu*self.mu*rp*rp-separation
        c1=(beta-k*d)**2-2*k*d+4*w*w*rp*rp-2*self.mu*self.mu*rp*d
        al=lambda n:mp.mpf((n+1)**2)
        be=lambda n:-2*n*(n-1)+b1*n+c0
        ga=lambda n:(n-1)*(n-2)+b2*(n-1)+c1
        ratio=mp.mpf(0);ratios=[]
        for n in range(self.terms,0,-1):
            ratio=-ga(n)/(be(n)+al(n)*ratio);ratios.append(ratio)
        return be(0)+al(0)*ratio,(rp,rm,w,k,beta,angular,list(reversed(ratios)))

    def radial(self,r,remove_exponential=False):
        with mp.workdps(self.dps):
            r=mp.mpf(float(r))
            if r<=self.rp:raise ValueError('Use an exterior radius')
            x=(r-self.rp)/(r-self.rm)
            series=mp.polyval(list(reversed(self.coefficients)),x)
            derivative=mp.polyval([n*self.coefficients[n] for n in range(len(self.coefficients)-1,0,-1)],x)
            factor=(r-self.rm)**self.beta
            if not remove_exponential:factor*=mp.exp(-self.k*r)
            value=factor*series
            slope=value*(self.beta/(r-self.rm)-self.k)+factor*derivative*(self.rp-self.rm)/(r-self.rm)**2
            # When remove_exponential=True the slope also has exp(+k r)
            # removed, but remains the derivative of the original R.
            return np.array([float(value),float(slope)])

    def angular(self,theta):
        x=np.cos(theta);st=np.sin(theta);S=np.zeros_like(x);Sp=np.zeros_like(x)
        for i,c in enumerate(self.angular_coefficients):
            ell=1+2*i;norm=np.sqrt((2*ell+1)/(4*np.pi)*np.exp(gammaln(ell)-gammaln(ell+2)))
            poly=lpmv(1,ell,x);lower=lpmv(1,ell-1,x) if ell>1 else np.zeros_like(x)
            S+=float(c)*norm*poly;Sp+=float(c)*norm*(ell*x*poly-(ell+1)*lower)/st
        return S,Sp

    def integrals(self,nradial=128,ntheta=48):
        # Laguerre integrates exp(-2k(r-r+)) analytically. No production
        # integration grid, cloud normalization, or radial ODE solver enters.
        x,weights=roots_laguerre(nradial);z,wt=np.polynomial.legendre.leggauss(ntheta)
        rp,a,w,k,mu=[float(v) for v in (self.rp,self.a,self.omega,self.k,self.mu)]
        r=rp+x/(2*k);R,Rp=np.array([self.radial(radius,True) for radius in r]).T
        S,Sp=self.angular(np.arccos(z));rr=r[:,None];sig=rr*rr+a*a*z*z;delta=(rr-rp)*(rr-float(self.rm))
        gtt=-((rr*rr+a*a)**2-a*a*delta*(1-z*z))/(sig*delta)
        gtp=-2*a*rr/(sig*delta);gpp=(delta-a*a*(1-z*z))/(sig*delta*(1-z*z))
        field2=R[:,None]**2*S*S
        energy=sig*((-gtt*w*w+gpp+mu*mu)*field2+delta/sig*Rp[:,None]**2*S*S+R[:,None]**2*Sp*Sp/sig)
        charge=2*sig*(-gtt*w+gtp)*field2
        prefactor=2*np.pi*np.exp(-2*k*rp)/(2*k)
        return prefactor*np.dot(weights,energy@wt),prefactor*np.dot(weights,charge@wt)
