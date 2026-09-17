"""Independent finite variation of the undensitized KG operator.
Real symmetric metric perturbations are deliberately not in Lorenz gauge.
No radial solver, reconstruction, or printed scalar-source formula is used.
"""
from pathlib import Path
import hashlib,json
import numpy as np
from environment_source import kerr_metric
OUT=Path(__file__).resolve().parents[1]/'docs/root_cause_followup_20260917'

def hmetric(r,t):
    z=np.zeros((4,4),dtype=np.result_type(r,t,float))
    z[0,0]=.07*np.sin(t)/r;z[1,1]=.03*r/(r+1)
    z[2,2]=.04*r*np.cos(t);z[3,3]=.06*r*np.sin(t)**2
    z[0,1]=z[1,0]=.025/(r+1)
    z[0,3]=z[3,0]=.013*np.sin(t)**2/r
    z[1,2]=z[2,1]=.011*np.sin(t);z[2,3]=z[3,2]=.017*np.sin(t)**2
    return z

def derivative(fun,r,t):
    d=np.zeros((4,4,4));d[1]=fun(r+1e-25j,t).imag/1e-25
    d[2]=fun(r,t+1e-25j).imag/1e-25
    return d

def geometry(g,dg):
    inv=np.linalg.inv(g);G=np.zeros((4,4,4))
    for i in range(4):
        for j in range(4):G[:,i,j]=.5*inv@(dg[i,:,j]+dg[j,:,i]-dg[:,i,j])
    return inv,G

def field(r,t,w,m):
    b=-.12+.03j;R=np.exp(b*r);S=np.sin(t)+.1*np.cos(t)**2
    Sp=np.cos(t)-.2*np.cos(t)*np.sin(t)
    Spp=-np.sin(t)-.2*(np.cos(t)**2-np.sin(t)**2)
    f=R*S;grad=np.array([-1j*w*f,b*f,R*Sp,1j*m*f])
    partial=np.empty((4,4),complex)
    partial[0]=-1j*w*grad;partial[:,0]=partial[0]
    partial[3]=1j*m*grad;partial[:,3]=partial[3]
    partial[1,1]=b*b*f;partial[2,2]=R*Spp
    partial[1,2]=partial[2,1]=b*R*Sp
    return f,grad,partial

def pair(z):return [float(z.real),float(z.imag)]
def main():
    rows=[]
    for a,r,t in [(0.,3.2,.7),(.88,3.2,.7),(.88,10.,1.1),(.88,45.,2.1)]:
        g=kerr_metric(r,t,a);dg=derivative(lambda r,t:kerr_metric(r,t,a),r,t)
        h=hmetric(r,t);dh=derivative(hmetric,r,t);inv,G=geometry(g,dg)
        di=np.array([-inv@d@inv for d in dg]);hu=inv@h@inv
        dhu=np.array([di[c]@h@inv+inv@dh[c]@inv+inv@h@di[c] for c in range(4)])
        dtrace=np.array([np.einsum('ij,ij->',di[c],h)+np.einsum('ij,ij->',inv,dh[c]) for c in range(4)])
        C=np.einsum('aab->b',dhu)+np.einsum('aad,db->b',G,hu)+np.einsum('bad,ad->b',G,hu)-.5*inv@dtrace
        f,grad,partial=field(r,t,.2963,1);H=partial-np.einsum('kij,k->ij',G,grad)
        target=-np.einsum('ij,ij->',hu,H)-C@grad
        base=np.einsum('ij,ij->',inv,H)-.3**2*f
        sequence=[]
        for eps in [1e-2,3e-3,1e-3,3e-4,1e-4,3e-5,1e-5]:
            def Q(lam):
                gi,Gamma=geometry(g+lam*h,dg+lam*dh)
                return np.einsum('ij,ij->',gi,partial-np.einsum('kij,k->ij',Gamma,grad))-.3**2*f
            fd=(Q(eps)-Q(-eps))/(2*eps)
            sequence.append(dict(epsilon=eps,finite_variation=pair(fd),relative_error=float(abs(fd-target)/abs(target))))
        best=min(x['relative_error'] for x in sequence)
        assert best<1e-8,(a,r,t,best)
        rows.append(dict(a=a,r=r,theta=t,analytic_variation=pair(target),lorenz_divergence=C.tolist(),
            background_KG=pair(base),Dyson_printed_extra=pair(-.5*np.einsum('ij,ij->',inv,h)*base),
            best_relative_error=best,finite_difference_sequence=sequence))
    result=dict(status='independent_metric_finite_variation_pass',
        identity='delta_g[(Box-mu^2)Phi] = -h^{ab}Hessian_ab - (nabla_a barh^{ab})partial_b Phi',
        limitation='Manufactured off-shell scalar and non-Lorenz h test the exact variational identity, not the actual particle metric.',
        rows=rows,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/'kg_metric_variation_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({str((r['a'],r['r'])):r['best_relative_error'] for r in rows}))
if __name__=='__main__':main()
