"""Independent coordinate linearized Ricci test of the retained Kerr metric.
Use only values and first/two partial derivatives; no reconstruction field
equation or tetrad curvature identity is substituted in the Ricci calculation.
"""
from pathlib import Path
import json,hashlib,concurrent.futures
import numpy as np
from lorenz_ghp import KerrGHP
from lorenz_jet import Jet
from lorenz_metric import nonstatic_metric
OUT=Path(__file__).resolve().parents[1]/'docs/root_cause_followup_20260917'
def val(z):return np.array([[x.value for x in row] for row in z])
def pd(g,f,k,harmonic=True):
    if harmonic:return g.partial(f,k)
    return f.derivative(k-1) if k in [1,2] else Jet(0,f.order)
def arrays(g,h):
    metric=val(g.g);inv=val(g.inv);hv=val(h)
    dg=np.array([val([[pd(g,z,k,False) for z in row] for row in g.g]) for k in range(4)])
    ddg=np.array([[val([[pd(g,pd(g,z,k,False),q,False) for z in row] for row in g.g]) for k in range(4)] for q in range(4)])
    dh=np.array([val([[pd(g,z,k) for z in row] for row in h]) for k in range(4)])
    ddh=np.array([[val([[pd(g,pd(g,z,k),q) for z in row] for row in h]) for k in range(4)] for q in range(4)])
    return metric,inv,hv,dg,ddg,dh,ddh

def audit(g,h):
    metric,inv,h,dg,ddg,dh,ddh=arrays(g,h)
    hu=inv@h@inv;di=np.array([-inv@d@inv for d in dg])
    dhu=np.array([di[q]@h@inv+inv@dh[q]@inv+inv@h@di[q] for q in range(4)])
    G=np.zeros((4,4,4),complex);dG=G.copy();ddG=np.zeros((4,4,4,4),complex)
    for b in range(4):
        for c in range(4):
            A=dg[b,:,c]+dg[c,:,b]-dg[:,b,c]
            B=dh[b,:,c]+dh[c,:,b]-dh[:,b,c]
            G[:,b,c]=.5*inv@A;dG[:,b,c]=.5*(-hu@A+inv@B)
            for q in range(4):
                dA=ddg[q,b,:,c]+ddg[q,c,:,b]-ddg[q,:,b,c]
                dB=ddh[q,b,:,c]+ddh[q,c,:,b]-ddh[q,:,b,c]
                ddG[q,:,b,c]=.5*(-dhu[q]@A-hu@dA+di[q]@B+inv@dB)
    terms=np.array([np.einsum('aabc->bc',ddG),-np.einsum('caba->bc',ddG),
        np.einsum('aad,dbc->bc',dG,G),np.einsum('aad,dbc->bc',G,dG),
        -np.einsum('acd,dba->bc',dG,G),-np.einsum('acd,dba->bc',G,dG)])
    ricci=terms.sum(axis=0)
    # Normalized null tetrad provides a reproducible physical component scale.
    E=np.array([[x.value for x in row] for row in g.tetrad]);proj=lambda x:E@x@E.T
    rr=proj(ricci);tt=np.array([proj(v) for v in terms]);hh=proj(h)
    norm=float(np.linalg.norm(rr));scale=float(sum(np.linalg.norm(v) for v in tt))
    dtrace=np.array([np.einsum('ij,ij->',di[q],h)+np.einsum('ij,ij->',inv,dh[q]) for q in range(4)])
    C=np.einsum('aab->b',dhu)+np.einsum('aad,db->b',G,hu)+np.einsum('bad,ad->b',G,hu)-.5*inv@dtrace
    return dict(ricci_tetrad_norm=norm,ricci_relative_to_term_norm=norm/(scale or 1),
        term_tetrad_norm_sum=scale,metric_tetrad_norm=float(np.linalg.norm(hh)),
        ricci_trace=[float(z) for z in [np.einsum('ij,ij->',inv,ricci).real,np.einsum('ij,ij->',inv,ricci).imag]],
        ricci_coordinate=[[[float(z.real),float(z.imag)] for z in row] for row in ricci],
        lorenz_coordinate_norm=float(np.linalg.norm(C)))

def gauge_check():
    g=KerrGHP(4.,.8,.88,omega=.015,m=1,order=6)
    xi=[Jet(0,6),g.r*g.theta.sin(),Jet(0,6),Jet(0,6)]
    h=[[sum(xi[k]*pd(g,g.g[i][j],k,False)+g.g[k][j]*g.partial(xi[k],i)+g.g[i][k]*g.partial(xi[k],j) for k in range(4)) for j in range(4)] for i in range(4)]
    result=audit(g,h)
    assert result['ricci_relative_to_term_norm']<1e-12,result
    return result

def worker(case):
    r,t,order=case;a=.8771530275949366;r0=42.1;pieces=[];summary=[]
    for ell in range(1,7):
        g,h=nonstatic_metric(r,t,r0,a,ell,1,order=order)
        pieces.append(h);summary.append(dict(ell=ell,**audit(g,h)))
    total=[[sum(h[i][j] for h in pieces) for j in range(4)] for i in range(4)]
    return dict(r=r,theta=t,order=order,orbit=r0,a=a,m_g=1,per_ell=summary,total=audit(g,total))
def main():
    fixture=gauge_check();rows=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:
        for row in pool.map(worker,[(3.,.71,10),(20.,1.1,10),(50.,.9,10),(20.,1.1,8)]):
            rows.append(row);print(row['r'],row['order'],row['total']['ricci_relative_to_term_norm'],flush=True)
            result=dict(status='local_vacuum_coordinate_linearized_Ricci_audit',gauge_fixture=fixture,rows=rows,
                limitations=['Off-orbit local vacuum equation only; it does not test delta-function point-particle normalization or global homogeneous gauge.',
                'Term-norm normalized residual measures local cancellation; absolute residual and h norm are also retained.'],
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
            (OUT/'metric_vacuum_ricci_audit.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
