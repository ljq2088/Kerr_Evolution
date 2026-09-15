"""All ten explicit 2023 nonstatic Lorenz tetrad components.

Unnormalized real-coordinate tetrad order: l+, l-, m+, m-.
Formula inputs are pure radial/angular potentials, not Kinnersley-weighted
fields. Derived from arXiv:2306.16459, spin-2/1/0 metric component equations.
No production tensor differentiation routine constructs these components.
"""
import numpy as np
from lorenz_jet import Jet
from paper_jump_basis import operators,radial_jet,angular_jet

PAIRS=((0,0),(1,1),(2,2),(3,3),(0,2),(0,3),(1,2),(1,3),(0,1),(2,3))
SPINS=(0,0,2,-2,1,-1,1,-1,0,0)


def scalar_tetrad(g,h,k):
    D,Dd,L,Ld=operators(g);rho=g.zeta.conjugate();rc=g.zeta
    r=g.r;a=g.a;w=g.omega;s=g.theta.sin();c=g.theta.cos();sig=g.sigma;de=g.delta
    Q=g.m/s-a*w*s;K=(r*r+a*a)*w-a*g.m
    mixed=-sig.derivative(0)*de*k.derivative(0)+sig.derivative(1)*k.derivative(1)
    return [-(D(r*D(h))/(1j*w)+4*D(D(k))),
       -(-Dd(r*Dd(h))/(1j*w)+4*Dd(Dd(k))),
       -Ld(-a*c/w*Ld(h)+4*Ld(k),-1),
       -L(a*c/w*L(h)+4*L(k),-1),
       -(sig/(2j*w)*D(Ld(h))+a/w*(r*s*D(h)+c*Ld(h))+4*rho*D(Ld(k))-4*(Ld(k)-1j*a*s*D(k)))/rho,
       -(sig/(2j*w)*D(L(h))-a/w*(r*s*D(h)+c*L(h))+4*rc*D(L(k))-4*(L(k)+1j*a*s*D(k)))/rc,
       -(-sig/(2j*w)*Dd(Ld(h))+a/w*(r*s*Dd(h)+c*Ld(h))+4*rc*Dd(Ld(k))-4*(Ld(k)+1j*a*s*Dd(k)))/rc,
       -(-sig/(2j*w)*Dd(L(h))-a/w*(r*s*Dd(h)+c*L(h))+4*rho*Dd(L(k))-4*(L(k)-1j*a*s*Dd(k)))/rho,
       (-(K/w-2*sig)*h-2*(D(de*Dd(k))+Dd(de*D(k)))-4*mixed/sig)/de,
       -(a*s*Q/w-2*sig)*h-2*(L(Ld(k),1)+Ld(L(k),1))+4*mixed/sig]


def vector_tetrad(g,ell,p0,p1):
    D,Dd,L,Ld=operators(g);rho=g.zeta.conjugate();rc=g.zeta;de=g.delta;sig=g.sigma
    Sm,lam=angular_jet(g,-1,ell);Sp,_=angular_jet(g,1,ell)
    Pm=radial_jet(g,-1,lam,p0,p1);p=(-1)**(ell+g.m)
    B=np.sqrt(lam*lam+4*g.a*g.m*g.omega-4*g.a*g.a*g.omega*g.omega)
    Pp=de*D(D(Pm))/B;calP=D(Pm)+p*Dd(Pp);calS=Ld(Sm,1)-p*L(Sp,1)
    ll=(calP*calS-sig.derivative(0)/sig*(Pm+p*Pp)*calS+p*calP*sig.derivative(1)/sig*(Sp-p*Sm))/de
    return [p*2/de*D(Pp,-1)*calS,2/de*Dd(Pm,-1)*calS,
        p*2*calP*Ld(Sp,-1),-2*calP*L(Sm,-1),
        p*(D(calP)-2/rho*calP)*Sp+p*Pp/de*(Ld(calS)-2*rho.derivative(1)/rho*calS),
        -(D(calP)-2/rc*calP)*Sm+p*Pp/de*(L(calS)-2*rc.derivative(1)/rc*calS),
        p*(Dd(calP)-2/rc*calP)*Sp+Pm/de*(Ld(calS)-2*rc.derivative(1)/rc*calS),
        -(Dd(calP)-2/rho*calP)*Sm+Pm/de*(L(calS)-2*rho.derivative(1)/rho*calS),ll,-de*ll]


def tensor_tetrad(g,ell,p0,p1):
    D,Dd,L,Ld=operators(g);rho=g.zeta.conjugate();rc=g.zeta;de=g.delta;sig=g.sigma
    Sm,lam=angular_jet(g,-2,ell);Sp,_=angular_jet(g,2,ell)
    Pm=radial_jet(g,-2,lam,p0,p1);aw=g.a*g.omega;p=(-1)**(ell+g.m)
    A=np.sqrt(lam**2*(lam+2)**2+8*aw*lam*((g.m-aw)*(5*lam+6)+12*aw)+144*aw**2*(g.m-aw)**2)
    Pp=de**2*D(D(D(D(Pm))))/(A+p*12j*g.omega)
    beta=1/(-6j*g.omega);psim=beta*Pm*Sm;psip=beta*Pp*Sp
    am=D(Ld(psim,2))/(3j*g.omega*rc);ap=Dd(L(psip,2))/(-3j*g.omega*rc)
    cm=D(D(Ld(Ld(psim,2),1)))/(24*g.omega**2)
    cp=Dd(Dd(L(L(psip,2),1)))/(24*g.omega**2);c=cm-cp
    def DR(f):return D(f)-2/rho*f
    def DdR(f):return Dd(f)-2/rho*f
    def LR(f):return L(f)-2*rho.derivative(1)/rho*f
    def LdR(f):return Ld(f)-2*rho.derivative(1)/rho*f
    def DC(f):return D(f)-2/rc*f
    def DdC(f):return Dd(f)-2/rc*f
    def LC(f):return L(f)-2*rc.derivative(1)/rc*f
    def LdC(f):return Ld(f)-2*rc.derivative(1)/rc*f
    hlpmp=-2*rc**2/de*(Dd(L(psip,2))-sig.derivative(0)/sig*L(psip,2)-sig.derivative(1)/sig*Dd(psip))
    hlpmp+=DR(Ld(c)+rc**2*Dd(ap))+LdR(D(c)-rc**2/de*L(ap,1))
    hlmmm=2*rc**2/de*(D(Ld(psim,2))-sig.derivative(0)/sig*Ld(psim,2)-sig.derivative(1)/sig*D(psim))
    hlmmm+=DdR(-rc**2*D(am)+L(c))+LR(rc**2/de*Ld(am,1)+Dd(c))
    hlpmm=DC(L(c)-rc**2*D(am))+LC(D(c)-rc**2/de*L(ap,1))
    hlmmp=DdC(Ld(c)+rc**2*Dd(ap))+LdC(Dd(c)+rc**2/de*Ld(am,1))
    ll=(sig*(D(de*Dd(c))+Dd(de*D(c)))-2*de*sig.derivative(0)*c.derivative(0)+2*sig.derivative(1)*c.derivative(1)
       +sig*(D(rc**2*Ld(am,1))-Dd(rc**2*L(ap,1)))-rc**2*sig.derivative(0)*(Ld(am,1)-L(ap,1))
       +rc**2*sig.derivative(1)*(Dd(ap)-D(am)))/(sig*de)
    return [-Pp*(p*Ld(Ld(Sm,2),1)+L(L(Sp,2),1))/(6*g.omega**2*de**2),
       -Pm*(Ld(Ld(Sm,2),1)+p*L(L(Sp,2),1))/(6*g.omega**2*de**2),
       -(p*D(D(Pm))+Dd(Dd(Pp)))*Sp/(6*g.omega**2),
       -(D(D(Pm))+p*Dd(Dd(Pp)))*Sm/(6*g.omega**2),
       hlpmp,hlpmm,hlmmp,hlmmm,ll,-de*ll]


def tetrad_vectors(r,t,a):
    de=r*r-2*r+a*a;s=np.sin(t)
    return np.array([[(r*r+a*a)/de,1,0,a/de],[-(r*r+a*a)/de,1,0,-a/de],
        [1j*a*s,0,1,1j/s],[-1j*a*s,0,1,-1j/s]],dtype=np.result_type(r,t,a,np.complex128))


def to_bl(components,r,t,a):
    tetrad=np.zeros((4,4),complex)
    for (i,j),c in zip(PAIRS,components):tetrad[i,j]=tetrad[j,i]=complex(c.value if hasattr(c,'value') else c)
    inverse=np.linalg.inv(tetrad_vectors(r,t,a))
    return inverse@tetrad@inverse.T


def project_metric(metric,r,t,a):
    V=tetrad_vectors(r,t,a);h=V@metric@V.T
    return np.array([h[i,j] for i,j in PAIRS])
