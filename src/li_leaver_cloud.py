"""Fixed-spin Leaver quasibound spectrum and finite-series radial profile.

Supports the two Li cloud states at a=.88, mu=.3. No mass normalization or
frequency freezing is hidden here. The output is an unnormalized complex mode.
"""
import mpmath as mp


def recurrence(a,mu,w,ell,m,terms,angular_terms=12):
    rp=1+mp.sqrt(1-a*a);rm=2-rp;gap=rp-rm
    q=-mp.sqrt(mu*mu-w*w)
    if mp.re(q)>=0:raise ValueError('Decaying branch required')
    beta=(mu*mu-2*w*w)/q-1;s=-1j*(2*rp*w-a*m)/gap
    low=abs(m)+(ell-abs(m))%2
    C=lambda l:mp.sqrt(mp.mpf(l*l-m*m)/(4*l*l-1)) if l>abs(m) else mp.mpf(0)
    mat=mp.matrix(angular_terms);c2=a*a*(w*w-mu*mu)
    for i in range(angular_terms):
        l=low+2*i;mat[i,i]=l*(l+1)-c2*(C(l+1)**2+C(l)**2)
        if i<angular_terms-1:mat[i,i+1]=mat[i+1,i]=-c2*C(l+1)*C(l+2)
    A=min(mp.eig(mat,left=False,right=False),key=lambda z:abs(z-ell*(ell+1)))
    sep=A+a*a*w*w-2*a*m*w
    def reduced(x):
        r=(rp-rm*x)/(1-x);delta=(r-rp)*(r-rm);K=(r*r+a*a)*w-a*m
        lp=s/x+beta/(1-x)+q*gap/(1-x)**2
        lpp=-s/x**2+beta/(1-x)**2+2*q*gap/(1-x)**3
        return (1-x)**2*(x*(lpp+lp*lp)+lp)+K*K/delta-mu*mu*r*r-sep
    x1,x2=mp.mpf('.25'),mp.mpf('.75');c1=(reduced(x2)-reduced(x1))/(x2-x1);c0=reduced(x1)-c1*x1
    defect=reduced(mp.mpf('.5'))-c0-c1/2
    b0=1+2*s;b1=-2*b0+2*(q*gap+beta);b2=b0-2*beta
    al=lambda n:(n+1)*(n+b0)
    be=lambda n:-2*n*(n-1)+b1*n+c0
    ga=lambda n:(n-1)*(n-2)+b2*(n-1)+c1
    ratios={};tail=mp.mpc(0)
    for n in range(terms,0,-1):
        tail=-ga(n)/(be(n)+al(n)*tail);ratios[n]=tail
    return dict(residual=be(0)+al(0)*tail,A=A,defect=defect,ratios=ratios,
        rp=rp,rm=rm,q=q,beta=beta,s=s)

class LiLeaverCloud:
    def __init__(self,ell=1,m=1,a='.88',mu='.3',terms=150,dps=64):
        if (ell,m) not in ((1,1),(2,2)) or str(a) not in ('.88','0.88') or str(mu) not in ('.3','0.3'):
            raise ValueError('Bounded to the two Li states at a=.88, mu=.3')
        if int(terms)!=terms or terms<20 or int(dps)!=dps or dps<40:raise ValueError('Invalid precision or truncation')
        self.dps,self.terms,self.ell,self.m=dps,int(terms),ell,m
        with mp.workdps(dps):
            self.a,self.mu=mp.mpf(str(a)),mp.mpf(str(mu))
            seed=mp.mpc('.296294','2e-9') if ell==1 else mp.mpc('.298451','1e-12')
            self.omega=mp.findroot(lambda w:recurrence(self.a,self.mu,w,ell,m,self.terms)['residual'],
                (seed,seed+mp.mpf('1e-7')),tol=mp.mpf(10)**(-dps+12),maxsteps=40)
            self.data=recurrence(self.a,self.mu,self.omega,ell,m,self.terms)
            self.coefficients=[mp.mpc(1)]
            for n in range(1,self.terms+1):self.coefficients.append(self.coefficients[-1]*self.data['ratios'][n])

    def radial_mp(self,r,derivative=0):
        with mp.workdps(self.dps):
            r=mp.mpf(str(r));d=self.data
            if r<=d['rp']:raise ValueError('Exterior radius required')
            def f(z):
                x=(z-d['rp'])/(z-d['rm'])
                # a0=1 with prefactor x^s exp(q*(r-r+))*((r-r-)/(r+-r-))^beta.
                return mp.exp(d['q']*(z-d['rp']))*((z-d['rm'])/(d['rp']-d['rm']))**d['beta']*x**d['s']*mp.polyval(list(reversed(self.coefficients)),x)
            return mp.diff(f,r,derivative)

    def radial_residual(self,r):
        with mp.workdps(self.dps):
            r=mp.mpf(str(r));d=self.data;R,Rp,Rpp=[self.radial_mp(r,k) for k in range(3)]
            delta=(r-d['rp'])*(r-d['rm']);K=(r*r+self.a*self.a)*self.omega-self.a*self.m
            V=K*K/delta-self.mu*self.mu*r*r-self.a*self.a*self.omega*self.omega+2*self.a*self.m*self.omega-d['A']
            terms=(delta*Rpp,2*(r-1)*Rp,V*R)
            return abs(sum(terms))/max(sum(abs(x) for x in terms),mp.mpf('1e-100'))
