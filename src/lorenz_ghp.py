"""Kinnersley tetrad and weighted GHP derivatives, BL (-+++) convention.

Spin coefficients are computed from covariant tetrad derivatives, using the
definitions in 2406.12510v3 Appendix A, rather than importing NP sign tables.
"""
import numpy as np
from copy import deepcopy
from functools import lru_cache
from lorenz_jet import Jet


@lru_cache(maxsize=64,typed=True)
def _geometry(r,theta,a,order,coefficient_dtype):
    # Geometry has no Fourier frequency or m dependence. The dtype is part
    # of the key because diagnostics may switch Taylor coefficient precision.
    template=KerrGHP.__new__(KerrGHP)
    template._initialize_geometry(r,theta,a,order)
    return template.__dict__


class KerrGHP:
    def __init__(self,r,theta,a,omega=0.,m=0,order=6):
        self.order,self.omega,self.m=order,omega,m
        # Jets are mutable. Never share cached arrays with a caller, and
        # preserve aliases within each independent geometry using deepcopy.
        self.__dict__.update(deepcopy(_geometry(r,theta,a,order,np.dtype(Jet.coefficient_dtype).str)))

    def _initialize_geometry(self,r,theta,a,order):
        r,t=Jet.variable(r,0,order),Jet.variable(theta,1,order)
        self.r,self.theta,self.a=r,t,a
        s,c=t.sin(),t.cos()
        self.zeta=r-1j*a*c
        sigma,delta=r*r+a*a*c*c,r*r-2*r+a*a
        self.sigma,self.delta=sigma,delta
        z=lambda:Jet(0,order)
        g=[[z() for _ in range(4)] for _ in range(4)]
        g[0][0]=-1+2*r/sigma
        g[0][3]=g[3][0]=-2*a*r*s*s/sigma
        g[1][1],g[2][2]=sigma/delta,sigma
        g[3][3]=(r*r+a*a+2*a*a*r*s*s/sigma)*s*s
        inv=[[z() for _ in range(4)] for _ in range(4)]
        inv[0][0]=-((r*r+a*a)**2-a*a*delta*s*s)/(sigma*delta)
        inv[0][3]=inv[3][0]=-2*a*r/(sigma*delta)
        inv[1][1],inv[2][2]=delta/sigma,1/sigma
        inv[3][3]=(delta-a*a*s*s)/(sigma*delta*s*s)
        self.g,self.inv=g,inv
        self.tetrad=[[(r*r+a*a)/delta,Jet(1,order),z(),a/delta],
                     [(r*r+a*a)/(2*sigma),-delta/(2*sigma),z(),a/(2*sigma)],
                     [1j*a*s/(2**.5*self.zeta.conjugate()),z(),
                      1/(2**.5*self.zeta.conjugate()),1j/(2**.5*s*self.zeta.conjugate())]]
        self.tetrad.append([v.conjugate() for v in self.tetrad[2]])
        self.cov=[[sum(g[i][j]*v[j] for j in range(4)) for i in range(4)] for v in self.tetrad]
        def geometric_derivative(f,index):
            return f.derivative(index-1) if index in (1,2) else z()
        self.gamma=[[[sum(inv[k][l]*(geometric_derivative(g[l][j],i)
                        +geometric_derivative(g[l][i],j)-geometric_derivative(g[i][j],l))/2
                        for l in range(4)) for j in range(4)] for i in range(4)] for k in range(4)]
        nabla=[[[geometric_derivative(v[j],i)-sum(self.gamma[k][i][j]*v[k] for k in range(4))
                 for j in range(4)] for i in range(4)] for v in self.cov]
        def contract(tensor,first,second):
            return sum(self.tetrad[first][i]*self.tetrad[second][j]*tensor[i][j]
                       for i in range(4) for j in range(4))
        self.sc={}
        for suffix,(l,n,mm,mb) in (('',(0,1,2,3)),('p',(1,0,3,2))):
            self.sc['beta'+suffix]=(contract(nabla[mm],mm,mb)-contract(nabla[l],mm,n))/2
            self.sc['eps'+suffix]=(contract(nabla[mm],l,mb)-contract(nabla[l],l,n))/2
            self.sc['rho'+suffix]=-contract(nabla[l],mb,mm)
            self.sc['tau'+suffix]=-contract(nabla[l],n,mm)
        self.psi2=-1/self.zeta**3

    def partial(self,f,index):
        if index==0:
            return -1j*self.omega*f
        if index==3:
            return 1j*self.m*f
        return f.derivative(index-1)

    def derivative(self,f,p,q,kind):
        """Return (value,p_new,q_new), automatically tracking GHP weights."""
        idx={'thorn':0,'thornp':1,'eth':2,'ethp':3}[kind]
        d=sum(v*self.partial(f,j) for j,v in enumerate(self.tetrad[idx]))
        sc=self.sc
        if idx==0:
            return d-p*sc['eps']*f-q*sc['eps'].conjugate()*f,p+1,q+1
        if idx==1:
            return d+p*sc['epsp']*f+q*sc['epsp'].conjugate()*f,p-1,q-1
        if idx==2:
            return d-p*sc['beta']*f+q*sc['betap'].conjugate()*f,p+1,q-1
        return d+p*sc['betap']*f-q*sc['beta'].conjugate()*f,p-1,q+1

    def scalar_wave(self,f):
        """Independent coordinate divergence form of Box for GHP weight (0,0)."""
        density=self.sigma*self.theta.sin()
        return sum(self.partial(density*sum(self.inv[i][j]*self.partial(f,j)
                                           for j in range(4)),i) for i in range(4))/density

    def teukolsky(self,f,spin):
        """Appendix B GHP operator, including scalar O=-Box/2 convention."""
        sc=self.sc
        p,q=2*spin,0
        if spin>=0:
            d1,d2,e1,e2='thorn','thornp','eth','ethp'
            rho,rhop,tau,taup=sc['rho'],sc['rhop'],sc['tau'],sc['taup']
            rhob,taubp=rho.conjugate(),taup.conjugate()
        else:
            d1,d2,e1,e2='thornp','thorn','ethp','eth'
            rho,rhop,tau,taup=sc['rhop'],sc['rho'],sc['taup'],sc['tau']
            rhob,taubp=rho.conjugate(),taup.conjugate()
        inner,pp,qq=self.derivative(f,p,q,d2)
        inner=inner-rhop*f
        first=self.derivative(inner,pp,qq,d1)[0]-(2*abs(spin)*rho+rhob)*inner
        inner,pp,qq=self.derivative(f,p,q,e2)
        inner=inner-taup*f
        second=self.derivative(inner,pp,qq,e1)[0]-(2*abs(spin)*tau+taubp)*inner
        curvature={0:1,1:0,2:3}[abs(spin)]
        return first-second-curvature*self.psi2*f
