"""Independent 2023 local sourced jump solve (nonstatic circular modes).

Development diagnostic: spin-2 curvature jumps use pybhpt; gauge jumps are
solved from point-particle metric matching, not the 2024 compact currents.
"""
from functools import lru_cache
import numpy as np
from scipy.linalg import lstsq
from pybhpt.swsh import Yslm
from environment_angular_diagnostic import DenseRealHarmonic
from environment_source import kerr_metric
from lorenz_ghp import KerrGHP
from lorenz_jet import Jet
from lorenz_metric import _homogeneous_radial_data
from lorenz_weyl import weyl_amplitudes
from paper_jump_basis import radial_jet,angular_jet,operators,scalar_components,projected_matrix


def orbit(r0,a):
    op=1/(r0**1.5+a);metric=kerr_metric(r0,np.pi/2,a)
    ut=1/np.sqrt(-metric[0,0]-2*op*metric[0,3]-op**2*metric[3,3])
    return op,ut,metric@np.array([ut,0.,0.,ut*op])


@lru_cache(maxsize=256)
def curvature_jumps(r0,a,m,ell):
    w=m/(r0**1.5+a);up,inn=weyl_amplitudes(r0,a,ell,m)[-2]
    lam,ru,du=_homogeneous_radial_data(-2,ell,m,a,w,r0,'Up')
    _,ri,di=_homogeneous_radial_data(-2,ell,m,a,w,r0,'In')
    # The paper's rescaled Weyl scalar is -4 times pybhpt's spin -2 field.
    # This convention bridge is checked against all published Weyl jumps.
    return lam,-4*(up*ru-inn*ri),-4*(up*du-inn*di)


def tensor_components(g,ell,j0,j1):
    D,Dd,L,Ld=operators(g);rhoc=g.zeta;rho=g.zeta.conjugate()
    Sm,lam=angular_jet(g,-2,ell);Sp,_=angular_jet(g,2,ell)
    Pm=radial_jet(g,-2,lam,j0,j1);aw=g.a*g.omega;p=(-1)**(ell+g.m)
    A=np.sqrt(lam**2*(lam+2)**2+8*aw*lam*((g.m-aw)*(5*lam+6)+12*aw)+144*aw**2*(g.m-aw)**2)
    Pp=g.delta**2*D(D(D(D(Pm))))/(A+p*12j*g.omega)
    beta=1/(-6j*g.omega);psim=beta*Pm*Sm;psip=beta*Pp*Sp
    am=D(Ld(psim,2))/(3j*g.omega*rhoc)
    ap=Dd(L(psip,2))/(-3j*g.omega*rhoc)
    cm=D(D(Ld(Ld(psim,2),1)))/(24*g.omega**2)
    cp=Dd(Dd(L(L(psip,2),1)))/(24*g.omega**2)
    c=cm-cp
    hll=-Pp*(p*Ld(Ld(Sm,2),1)+L(L(Sp,2),1))/(6*g.omega**2*g.delta**2)
    hmm=-(p*D(D(Pm))+Dd(Dd(Pp)))*Sp/(6*g.omega**2)
    v=Ld(c)+rhoc**2*Dd(ap)
    u=D(c)-rhoc**2/g.delta*L(ap,1)
    hlm=-2*rhoc**2/g.delta*(Dd(L(psip,2))-g.sigma.derivative(0)/g.sigma*L(psip,2)-g.sigma.derivative(1)/g.sigma*Dd(psip))
    hlm+=D(v)-2/rho*v+Ld(u)-2*rho.derivative(1)/rho*u
    return [hll,hmm,rho*hlm]


@lru_cache(maxsize=64)
def trace_angular_data(r0,a,m,L):
    x,w=np.polynomial.legendre.leggauss(max(32,2*L+8));theta=np.arccos(x)
    harmonics=[DenseRealHarmonic(0,l,m,a*m/(r0**1.5+a)) for l in range(abs(m),L+1)]
    S=np.array([h(theta) for h in harmonics]);gamma=2*np.pi*np.einsum('ik,jk,k->ij',S,S,w*x*x)
    gamma[abs(gamma)<1e-14]=0.
    return np.array([h.eigenvalue for h in harmonics]),np.array([h(np.pi/2) for h in harmonics]),gamma


def known_components(g,r0,L):
    _,ut,_=orbit(r0,g.a);n=g.order
    lam,eq,gamma=trace_angular_data(r0,g.a,g.m,L)
    sj=[angular_jet(g,0,l)[0] for l in range(abs(g.m),L+1)]
    h=Jet(0.,n);k=Jet(0.,n);tensor=[Jet(0.,n) for _ in range(3)]
    for i,ell in enumerate(range(abs(g.m),L+1)):
        hr=radial_jet(g,0,lam[i],0.,-16*np.pi*eq[i]/(ut*float(g.delta.value.real)))
        h+=hr*sj[i]
        # Local diagonal particular solution with zero value/first derivative.
        kr=Jet(0.,n);K=(g.r*g.r+g.a*g.a)*g.omega-g.a*g.m
        V=K*K/g.delta-lam[i];src=(g.r*g.r+g.a*g.a*gamma[i,i])*hr/2
        for j in range(n-1):
            rhs=(src-2*(g.r-1)*kr.derivative(0)-V*kr)/g.delta
            kr.c[j+2,0]=rhs.c[j,0]/((j+1)*(j+2))
        k+=kr*sj[i]
        for j in range(len(sj)):
            if i!=j and gamma[i,j]!=0:
                k+=g.a*g.a*gamma[i,j]/(2*(lam[i]-lam[j]))*hr*sj[j]
        if ell>=2:
            _,j0,j1=curvature_jumps(r0,g.a,g.m,ell)
            term=tensor_components(g,ell,j0,j1)
            tensor=[v+term[j] for j,v in enumerate(tensor)]
    scalar=scalar_components(g,h,k)
    return [tensor[j]+scalar[j] for j in range(3)]


def solve_jumps(r0,a,m,L,q=20,order=10):
    A,degrees,labels=projected_matrix(r0,a,m,L,L,q)
    x,w=np.polynomial.legendre.leggauss(q);theta=np.arccos(x);spins=(0,2,1)
    angular=np.array([[Yslm(s,j,m,theta) if j>=max(abs(m),abs(s)) else np.zeros_like(theta) for s in spins] for j in degrees])
    samples=[]
    for t in theta:
        g=KerrGHP(r0,t,a,omega=m/(r0**1.5+a),m=m,order=order)
        h=known_components(g,r0,L)
        samples.append([[v.value for v in h],[v.derivative(0).value for v in h]])
    known=2*np.pi*np.einsum('jfk,kdf,k->djf',angular,np.array(samples),w)
    op,ut,ucov=orbit(r0,a);delta=r0*r0-2*r0+a*a
    lp=np.array([(r0*r0+a*a)/delta,1.,0.,a/delta])
    mp=np.array([1j*a,0.,1.,1j])
    targets=np.array([(ucov@lp)**2,(ucov@mp)**2,r0*(ucov@lp)*(ucov@mp)])
    target=np.zeros_like(known)
    target[1]=np.array([[-16*np.pi*targets[c]*Yslm(s,j,m,np.pi/2)/(ut*delta) if j>=max(abs(m),abs(s)) else 0 for c,s in enumerate(spins)] for j in degrees])
    mask=np.zeros(A.shape[:-1],bool);mask[:,:,0]=mask[:,:,1]=True
    if abs(m)==1:mask[:,0,1]=False;mask[:,0,2]=True
    mat=A[mask];b=(target-known)[mask]
    rowscale=np.linalg.norm(mat,axis=1);mat=mat/rowscale[:,None];b=b/rowscale
    colscale=np.linalg.norm(mat,axis=0);mat=mat/colscale[None,:]
    y,_,rank,_=lstsq(mat,b,lapack_driver='gelsy');solution=y/colscale
    if rank!=len(labels):raise RuntimeError('Matching system rank deficient')
    residual=np.einsum('djfc,c->djf',A,solution)+known-target
    return dict(r0=r0,a=a,m=m,L=L,q=q,order=order,labels=labels,solution=solution,
      condition_number=float(np.linalg.cond(np.asarray(mat,complex))),residual=residual,target=target,
      selected_residual=float(np.max(abs(residual[mask]))),
      unused_residual=float(np.max(abs(residual[~mask]))))
