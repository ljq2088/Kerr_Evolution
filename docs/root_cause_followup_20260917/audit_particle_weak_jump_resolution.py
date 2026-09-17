"""Independent distributional Einstein shell coefficients of a truncated metric.
Retains [h] delta-prime and all associated trace/weight derivatives. No production
source formula is used as reference. Diagnostic only; no production changes.
"""
from pathlib import Path
import sys,json,time,hashlib,concurrent.futures
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from environment_source import kerr_metric,connection
from lorenz_ghp import KerrGHP
from lorenz_metric import nonstatic_metric
from source_provenance import source_fingerprint
BASE=Path(__file__).resolve().parent
A=.8771530275949366;RP=42.1;OMG=1/(RP**1.5+A)
OFFSETS=[.001,.0005]

def geometry(r,t,a):
    g=kerr_metric(r,t,a);inv,G=connection(r,t,a);dg=np.zeros((4,4,4),complex)
    dg[1]=kerr_metric(r+1e-25j,t,a).imag/1e-25;dg[2]=kerr_metric(r,t+1e-25j,a).imag/1e-25
    di=np.array([-inv@q@inv for q in dg]);return g,inv,G,dg,di

def delta_connection(inv,dg,h,dh):
    hu=inv@h@inv;out=np.zeros((4,4,4),complex)
    for b in range(4):
        for c in range(4):
            av=dg[b,:,c]+dg[c,:,b]-dg[:,b,c]
            bv=dh[b,:,c]+dh[c,:,b]-dh[:,b,c]
            out[:,b,c]=.5*(-hu@av+inv@bv)
    return out

def shell(r,t,a,left,right):
    g,inv,G,dg,di=geometry(r,t,a);hl,dhl=left;hr,dhr=right
    J=hr-hl;dJ=dhr-dhl;n=np.array([0.,1.,0.,0.])
    def coeff(j,iv):
        out=np.zeros((4,4,4),complex)
        for b in range(4):
            for c in range(4):out[:,b,c]=.5*iv@(n[b]*j[:,c]+n[c]*j[:,b]-n*j[b,c])
        return out
    aa=coeff(J,inv);da=np.zeros((4,4,4,4),complex)
    # A is evaluated on the shell and extended constantly along r before
    # differentiating its delta distribution. Only tangential derivatives act.
    for q in [0,2,3]:da[q]=coeff(dJ[q],inv)+coeff(J,di[q])
    jump=delta_connection(inv,dg,hr,dhr)-delta_connection(inv,dg,hl,dhl)
    C=jump[1]-np.outer(np.einsum('aba->b',jump),n)
    C+=np.einsum('aabc->bc',da)-np.einsum('caba->bc',da)
    C+=np.einsum('aad,dbc->bc',aa,G)+np.einsum('aad,dbc->bc',G,aa)-np.einsum('acd,dba->bc',aa,G)-np.einsum('acd,dba->bc',G,aa)
    D=aa[1]-np.outer(np.einsum('aba->b',aa),n)
    tr=lambda x:np.einsum('ij,ij->',inv,x)
    P_D=D-.5*g*tr(D)
    # P(r) D delta' = P0 D delta' - P'_0 D delta.
    P_C=C-.5*g*tr(C)+.5*dg[1]*tr(D)+.5*g*np.einsum('ij,ij->',di[1],D)
    sig=r*r+a*a*np.cos(t)**2
    # A radial weak test constant near the shell has (Sigma*P_D)' term.
    weighted_C=sig*P_C-2*r*P_D;weighted_D=sig*P_D
    # Lorenz delta coefficient, independently exposing metric-jump constraints.
    lu=inv@J@inv;L=lu[1]-.5*inv[1]*tr(J)
    return dict(weighted_delta=weighted_C,weighted_delta_prime=weighted_D,jump_h=J,
        jump_dh=dJ,ricci_delta=C,ricci_delta_prime=D,einstein_delta=P_C,einstein_delta_prime=P_D,
        lorenz_delta=L)

def fixture():
    r,t,a=4.,.8,.88;g,inv,G,dg,di=geometry(r,t,a)
    z=np.zeros_like(g);dz=np.zeros_like(dg)
    out=shell(r,t,a,(z,dz),(g,dg))
    # h_ab = g_ab Theta(r-r0): delta G=-nabla_a nabla_b Theta+g_ab Box Theta.
    n=np.array([0.,1.,0.,0.]);traceG=np.einsum('ab,ab->',inv,G[1])
    exactD=-np.outer(n,n)+g*inv[1,1]
    exactC=G[1]-g*traceG-dg[1]*inv[1,1]-g*di[1,1,1]
    errors=dict(delta=float(np.linalg.norm(out['einstein_delta']-exactC)/np.linalg.norm(exactC)),
        delta_prime=float(np.linalg.norm(out['einstein_delta_prime']-exactD)/np.linalg.norm(exactD)))
    assert max(errors.values())<1e-12,errors
    return errors

def jets(g,h,dist):
    val=lambda q:np.array([[x.value for x in row] for row in q])
    dh=[[ [g.partial(h[i][j],k) for j in range(4)] for i in range(4)] for k in range(4)]
    drdh=[[[g.partial(dh[k][i][j],1) for j in range(4)] for i in range(4)] for k in range(4)]
    hv=val(h);d=np.array([val(q) for q in dh]);dd=np.array([val(q) for q in drdh])
    # Values quadratic / all first derivatives linear extrapolated to r0.
    h0=hv-dist*d[1]+.5*dist*dist*dd[1];d0=d-dist*dd
    return h0,d0

def worker(task):
    theta,ellmax=task
    rows=[]
    for eps in OFFSETS:
        sides=[]
        for sign in [-1,1]:
            pieces=[]
            for ell in range(1,ellmax+1):
                g,h=nonstatic_metric(RP+sign*eps,float(theta),RP,A,ell,1,order=8);pieces.append(h)
            total=[[sum(h[i][j] for h in pieces) for j in range(4)] for i in range(4)]
            sides.append(jets(g,total,sign*eps))
        rows.append(dict(offset=eps,**shell(RP,float(theta),A,*sides)))
    return float(theta),rows

def encode(obj):
    if isinstance(obj,dict):return {k:encode(v) for k,v in obj.items()}
    if isinstance(obj,(tuple,list)):return [encode(v) for v in obj]
    if isinstance(obj,np.ndarray):return encode(obj.tolist())
    if isinstance(obj,complex):return [obj.real,obj.imag]
    if isinstance(obj,np.generic):return obj.item()
    return obj

def main():
    import argparse
    pa=argparse.ArgumentParser();pa.add_argument('--L',type=int,required=True);pa.add_argument('--q',type=int,required=True);args=pa.parse_args()
    fp=source_fingerprint();unit=fixture();x,w=np.polynomial.legendre.leggauss(args.q);theta=np.arccos(x)
    result=dict(status='running',parameters=dict(a=A,rp=RP,mg=1,Lmetric=args.L,angular_order=args.q,jet_order=8,offsets=OFFSETS,workers=4),fixture=unit,source_provenance=fp,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),samples=[])
    out=BASE/f'particle_weak_jump_L{args.L}_q{args.q}.json'
    save=lambda:out.write_text(json.dumps(encode(result),indent=2,allow_nan=False)+'\n')
    started=time.perf_counter();save();values={}
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
        for t,rows in pool.map(worker,[(t,args.L) for t in theta]):
            values[t]=rows;result['samples'].append(dict(theta=t,rows=rows));save();print('ANGLE_COMPLETE',len(values),args.q,'elapsed',time.perf_counter()-started,flush=True)
    g=kerr_metric(RP,np.pi/2,A);ut=1/np.sqrt(-g[0,0]-2*OMG*g[0,3]-OMG**2*g[3,3]);uc=g@np.array([ut,0.,0.,OMG*ut]);particle=4*np.outer(uc,uc)/ut
    projections=[]
    # Even low-order Legendre test functions all have nonzero equatorial value;
    # odd tests provide zero-source parity diagnostics. No complex conjugate of h.
    for ie,eps in enumerate(OFFSETS):
        C=np.array([values[float(t)][ie]['weighted_delta'] for t in theta]);D=np.array([values[float(t)][ie]['weighted_delta_prime'] for t in theta]);Jh=np.array([values[float(t)][ie]['jump_h'] for t in theta])
        for k in range(5):
            coeff=np.zeros(k+1);coeff[-1]=1;test=np.polynomial.legendre.legval(x,coeff);test0=np.polynomial.legendre.legval(0.,coeff)
            actual=np.einsum('t,tbc->bc',w*test,C);dprime=np.einsum('t,tbc->bc',w*test,D);expected=particle*test0
            projections.append(dict(offset=eps,test_legendre_degree=k,test_at_particle=float(test0),actual_delta=actual,expected_delta=expected,delta_prime=dprime,
                relative_error_to_particle_norm=float(np.linalg.norm(actual-expected)/np.linalg.norm(particle)),
                delta_prime_over_M_particle_norm=float(np.linalg.norm(dprime)/np.linalg.norm(particle)),
                symmetry_error=float(np.linalg.norm(actual-actual.T)/np.linalg.norm(particle))))
    result.update(status='finite_L_distributional_shell_audit_complete',elapsed_seconds=time.perf_counter()-started,particle_u_covariant=uc,particle_ut=ut,particle_reference=particle,projections=projections,
        distribution_convention='Linearized Einstein tensor integrated against smooth radial test f constant near r0 and polar Pk(cos theta), with Sigma measure. Delta-prime tested separately by f derivative.',
        expected='8pi integral Sigma T_m dr dx = 4 u_a u_b/u^t * Pk(0), q=1 and Fourier coefficient 1/(2pi).',
        limitations=['Finite L and q as recorded; compare multiple outputs before interpreting convergence.',
            'One-sided fields extrapolated from two finite offsets using value/first/second jets; derivative extrapolation is first order.',
            'All h jump induced delta-prime, tangential-derivative, trace-projection and Sigma weight derivative terms retained.',
            'This is a geometric weak distribution check, not an application of reconstruction source formulas.'])
    assert source_fingerprint()==fp
    save();print('COMPLETE',[(p['offset'],p['test_legendre_degree'],p['relative_error_to_particle_norm'],p['delta_prime_over_M_particle_norm']) for p in projections],flush=True)
if __name__=='__main__':main()
