"""2023 sourced matching using all ten explicit tensor component formulas.

Four equations per degree determine gauge jumps; the remaining tensor
conditions are held out. No comparison with environmental flux enters.
"""
import numpy as np
from scipy.linalg import lstsq
from pybhpt.swsh import Yslm
from lorenz_ghp import KerrGHP
from lorenz_jet import Jet
from paper_jump_basis import radial_jet,angular_jet
from paper_full_tetrad import scalar_tetrad,vector_tetrad,tensor_tetrad,SPINS,tetrad_vectors,PAIRS
from paper_sourced_matching import orbit,trace_angular_data,curvature_jumps
from environment_source import kerr_metric


def weighted(g,h):
    rho=g.zeta.conjugate();rc=g.zeta
    return [*h[:4],rho*h[4],rc*h[5],rc*h[6],rho*h[7],g.sigma*g.delta*h[8],(g.delta*h[8]+h[9])/g.sigma]


def known(g,r0,L):
    _,ut,_=orbit(r0,g.a);n=g.order;lam,eq,gamma=trace_angular_data(r0,g.a,g.m,L)
    angular=[angular_jet(g,0,ell)[0] for ell in range(abs(g.m),L+1)]
    h=Jet(0.,n);k=Jet(0.,n);tensor=[Jet(0.,n) for _ in range(10)]
    for i,ell in enumerate(range(abs(g.m),L+1)):
        hr=radial_jet(g,0,lam[i],0.,-16*getattr(g,'normalization_pi',np.pi)*eq[i]/(ut*g.delta.value.real));h+=hr*angular[i]
        kr=Jet(0.,n);K=(g.r*g.r+g.a*g.a)*g.omega-g.a*g.m
        V=K*K/g.delta-lam[i];src=(g.r*g.r+g.a*g.a*gamma[i,i])*hr/2
        for j in range(n-1):
            rhs=(src-2*(g.r-1)*kr.derivative(0)-V*kr)/g.delta
            kr.c[j+2,0]=rhs.c[j,0]/((j+1)*(j+2))
        k+=kr*angular[i]
        for j in range(len(angular)):
            if j!=i and gamma[i,j]!=0:k+=g.a*g.a*gamma[i,j]/(2*(lam[i]-lam[j]))*hr*angular[j]
        if ell>=2:
            _,j0,j1=curvature_jumps(r0,g.a,g.m,ell)
            term=tensor_tetrad(g,ell,j0,j1)
            tensor=[v+term[j] for j,v in enumerate(tensor)]
    sc=scalar_tetrad(g,h,k)
    return weighted(g,[tensor[j]+sc[j] for j in range(10)])


def solve(r0,a,m,L,q=24,order=10,extended_solve=False):
    if m<1:raise ValueError('Diagnostic uses positive mg and obtains the negative branch by conjugation')
    pi=np.pi;yslm=Yslm
    if extended_solve:
        from paper_precise_angular import precise_grid,precise_spherical
        r0,a=np.longdouble(r0),np.longdouble(a)
        theta,w,pi=precise_grid(q);yslm=precise_spherical
    else:
        x,w=np.polynomial.legendre.leggauss(q);theta=np.arccos(x)
    degrees=list(range(m,L+1))
    angular=np.array([[yslm(s,j,m,theta) if j>=max(abs(m),abs(s)) else np.zeros_like(theta) for s in SPINS] for j in degrees])
    labels=[(ell,kind,datum) for ell in degrees for kind in ('spin1','kappa') for datum in (0,1)]
    values=[];constants=[]
    for n,t in enumerate(theta):
        g=KerrGHP(r0,t,a,omega=m/(r0**1.5+a),m=m,order=order);columns=[];g.normalization_pi=pi
        for ell,kind,datum in labels:
            if kind=='spin1':h=vector_tetrad(g,ell,float(datum==0),float(datum==1))
            else:
                S,lam=angular_jet(g,0,ell);k=radial_jet(g,0,lam,float(datum==0),float(datum==1))*S
                h=scalar_tetrad(g,Jet(0.,order),k)
            h=weighted(g,h)
            columns.append([[v.value for v in h],[v.derivative(0).value for v in h]])
        values.append(columns);h=known(g,r0,L)
        constants.append([[v.value for v in h],[v.derivative(0).value for v in h]])
        print('matching angular node',n+1,'/',q,flush=True)
    A=2*pi*np.einsum('jfk,kcdf,k->djfc',angular,np.array(values),w)
    known_proj=2*pi*np.einsum('jfk,kdf,k->djf',angular,np.array(constants),w)
    _,ut,ucov=orbit(r0,a);g0=kerr_metric(r0,pi/2,a);de=r0*r0-2*r0+a*a
    hsource=-16*pi*(np.outer(ucov,ucov)+g0/2)/(ut*de)
    V=tetrad_vectors(r0,pi/2,a);H=V@hsource@V.T
    source=np.array([H[i,j] for i,j in PAIRS]);source[4:8]*=r0;source[8]*=r0*r0*de
    source[9]=-16*pi/(ut*de)  # u.u + 4/2 = 1, no matrix inverse needed
    target=np.zeros_like(known_proj)
    target[1]=np.array([[source[c]*yslm(s,j,m,pi/2) if j>=max(abs(m),abs(s)) else 0 for c,s in enumerate(SPINS)] for j in degrees])
    mask=np.zeros(A.shape[:-1],bool);mask[:,:,0]=mask[:,:,2]=True
    if m==1:mask[:,0,2]=False;mask[:,0,4]=True
    B=A[mask];b=(target-known_proj)[mask];rs=np.linalg.norm(B,axis=1);B=B/rs[:,None];b=b/rs
    cs=np.linalg.norm(B,axis=0);B=B/cs[None,:];solution,_,rank,_=lstsq(B,b,lapack_driver='gelsy');solution=solution/cs
    if rank!=len(labels):raise RuntimeError('Rank-deficient physical matching')
    # scipy LAPACK downcasts complex256 to complex128. Keep that result as
    # a control and optionally solve the same scaled system at 50 digits.
    def errors(solution):
        residual=np.einsum('djfc,c->djf',A,solution)+known_proj-target
        scale=np.maximum(np.max(abs(target[1]),axis=0),1e-10)
        rows=[dict(j=j,absolute_max=float(np.max(abs(residual[:,i]))),
            scaled_max=float(np.max(abs(residual[:,i])/scale))) for i,j in enumerate(degrees)]
        return residual,rows
    qr_solution=solution.copy();qr_residual,qr_rows=errors(solution)
    if extended_solve:
        import mpmath as mp
        def convert(z):
            fmt=lambda x:np.format_float_scientific(np.longdouble(x),precision=35,unique=False)
            return mp.mpc(fmt(np.real(z)),fmt(np.imag(z)))
        with mp.workdps(50):
            precise=mp.lu_solve(mp.matrix([[convert(z) for z in row] for row in B]),mp.matrix([convert(z) for z in b]))
            solution=np.array([np.clongdouble(np.longdouble(str(z.real)))+1j*np.longdouble(str(z.imag)) for z in precise])/cs
    residual,rows=errors(solution)
    return dict(r0=float(r0),a=float(a),m=m,L=L,q=q,order=order,labels=labels,solution=solution,
        residual=residual,target=target,selected_mask=mask,condition_number=float(np.linalg.cond(np.asarray(B,complex))),
        degree_errors=rows,extended_solve=extended_solve,qr_solution=qr_solution,qr_degree_errors=qr_rows,selected_residual=float(np.max(abs(residual[mask]))),
        unused_residual=float(np.max(abs(residual[~mask]))))
