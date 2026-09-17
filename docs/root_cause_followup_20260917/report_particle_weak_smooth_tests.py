"""Smooth sphere weak tests from saved one-sided metric shells, no metric rerun."""
from pathlib import Path
import json,hashlib
import numpy as np
B=Path(__file__).resolve().parent

def matrix(z):
    a=np.asarray(z)
    return a[...,0]+1j*a[...,1] if a.shape==(4,4,2) else a.astype(complex)
def enc(z):
    if isinstance(z,dict):return {k:enc(v) for k,v in z.items()}
    if isinstance(z,list):return [enc(v) for v in z]
    if isinstance(z,np.ndarray):return enc(z.tolist())
    if isinstance(z,complex):return [z.real,z.imag]
    if isinstance(z,np.generic):return z.item()
    return z

def main():
    runs=[];inputs={}
    for path in [B/'particle_weak_jump_audit.json']+sorted(B.glob('particle_weak_jump_L*_q*.json')):
        data=json.loads(path.read_text())
        if data['status']!='finite_L_distributional_shell_audit_complete':continue
        inputs[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        q=data['parameters']['angular_order'];x,w=np.polynomial.legendre.leggauss(q);theta=np.arccos(x);st=np.sin(theta)
        assert np.max(abs(np.array([r['theta'] for r in data['samples']])-theta))<1e-13
        target=matrix(data['particle_reference']);scale=np.linalg.norm(target)
        rows=[]
        for ie,offset in enumerate(data['parameters']['offsets']):
            C=np.array([matrix(v['rows'][ie]['weighted_delta']) for v in data['samples']]);D=np.array([matrix(v['rows'][ie]['weighted_delta_prime']) for v in data['samples']])
            # Smooth vector basis (dt, dr, sin(theta)dtheta, dphi). Unlike dtheta,
            # sin(theta)dtheta=-grad_{unit sphere}(cos theta) is globally smooth.
            V=np.ones((q,4));V[:,2]=st
            C=C*V[:,:,None]*V[:,None,:];D=D*V[:,:,None]*V[:,None,:]
            for k in range(5):
                coeff=np.zeros(k+1);coeff[-1]=1
                test=st*np.polynomial.legendre.legval(x,coeff);at0=np.polynomial.legendre.legval(0.,coeff)
                actual=np.einsum('t,tbc->bc',w*test,C);dprime=np.einsum('t,tbc->bc',w*test,D);expected=target*at0
                components={name:dict(actual=actual[i,j],expected=expected[i,j],relative_error=abs(actual[i,j]-expected[i,j])/abs(target[i,j]),delta_prime_over_M_reference=abs(dprime[i,j])/abs(target[i,j])) for name,(i,j) in {'tt':(0,0),'tphi':(0,3),'phiphi':(3,3)}.items()}
                rows.append(dict(offset=offset,test_polynomial_degree=k,test_at_particle=float(at0),actual_delta=actual,expected_delta=expected,delta_prime=dprime,
                    relative_error_to_particle_norm=float(np.linalg.norm(actual-expected)/scale),delta_prime_over_M_particle_norm=float(np.linalg.norm(dprime)/scale),components=components))
        runs.append(dict(file=path.name,parameters=data['parameters'],rows=rows))
        print(path.name,[(v['test_polynomial_degree'],v['relative_error_to_particle_norm'],v['delta_prime_over_M_particle_norm']) for v in rows if v['offset']==.0005],flush=True)
    result=dict(status='smooth_weak_test_projection_complete',test='sin(theta) Pk(cos theta) exp(-i phi); tensor tests built from smooth (dt,dr,sin(theta)dtheta,dphi) vector basis',
        explanation='sin(theta) exp(-i phi) is a smooth ell=1 scalar; multiplying by a polynomial preserves smoothness. For degree0/2/4 these span scalar ell<=1/3/5 at m1. Delta-prime terms retained.',
        normalization='8pi integral Sigma T_m dr dx =4 u_a u_b/u^t times Pk(0), because sin(pi/2)=1 and all reference theta components vanish.',
        original_Pk_test_warning='Unweighted Pk(cos theta) exp(-i phi) is not smooth at sphere axes; original Pk diagnostics cannot alone establish particle source mismatch.',
        inputs_sha256=inputs,runs=runs,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (B/'particle_weak_smooth_tests.json').write_text(json.dumps(enc(result),indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
