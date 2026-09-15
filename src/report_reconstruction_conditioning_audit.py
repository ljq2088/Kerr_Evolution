"""Local term-cancellation audit, actual Kerr and Schwarzschild rp=20 mg=1.

Splits scalar CKY/kappa/compact-chi and tensor AAB/DKW-vector/chi. Norms are
measured after transformation to a local orthonormal ZAMO frame. These are
terminal-sum condition estimates, not a rigorous bound for all internal radial
or Taylor operations and not a source projection or final flux validation.
"""
import json,hashlib,time
from pathlib import Path
import numpy as np
from lorenz_ghp import KerrGHP
from lorenz_metric import spin1_metric,spin2_metric,homogeneous_field_jet
from lorenz_kappa import trace_field_jet,kappa_jet
from lorenz_spin1 import cky_tensor
from lorenz_chi import chi_amplitudes
from lorenz_tensor import vector_covariant_derivative
from environment_source import kerr_metric

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'docs/environment_reproduction/reconstruction_conditioning_audit.json'

def values(matrix):return np.array([[v.value for v in row] for row in matrix],dtype=complex)

def gauge_piece(g,xi):
    dx=vector_covariant_derivative(g,xi)
    return values([[-(dx[i][j]+dx[j][i])*1j/g.omega for j in range(4)] for i in range(4)])

def pieces(r,theta,r0,a,ell,m,order=6):
    g=KerrGHP(r,theta,a,omega=m/(r0**1.5+a),m=m,order=order)
    trace=trace_field_jet(g,r0,ell);kappa=kappa_jet(g,r0,ell)
    bc,index=('In',1) if r<r0 else ('Up',0)
    chi=homogeneous_field_jet(g,0,ell,chi_amplitudes(r0,a,ell,m)[index],bc)
    ky=cky_tensor(g)
    cky=[sum(ky[i][j]*g.inv[j][k]*g.partial(trace,k)/2 for j in range(4) for k in range(4)) for i in range(4)]
    out={'s0_CKY':gauge_piece(g,cky),
         's0_kappa':gauge_piece(g,[g.partial(kappa,i) for i in range(4)]),
         's0_compact_chi':gauge_piece(g,[-g.partial(chi,i) for i in range(4)])}
    g,xi1=spin1_metric(r,theta,r0,a,ell,m,order=order,return_vector=True,full_current=True)
    out['s1']=gauge_piece(g,xi1)
    if ell>=2:
        g,aab,vector,chi=spin2_metric(r,theta,r0,a,ell,m,order=order,return_parts=True)
        out.update(s2_AAB=values(aab)*1j/g.omega,s2_DKW_vector=gauge_piece(g,vector),
                   s2_chi=gauge_piece(g,[-g.partial(chi,i) for i in range(4)]))
    return out

def frame(r,theta,a):
    sigma=r*r+a*a*np.cos(theta)**2;delta=r*r-2*r+a*a
    area=(r*r+a*a)**2-a*a*delta*np.sin(theta)**2
    lapse=np.sqrt(sigma*delta/area);drag=2*a*r/area
    e=np.zeros((4,4));e[0,0]=1/lapse;e[0,3]=drag/lapse
    e[1,1]=np.sqrt(delta/sigma);e[2,2]=1/np.sqrt(sigma)
    e[3,3]=np.sqrt(sigma/area)/np.sin(theta)
    np.testing.assert_allclose(e@kerr_metric(r,theta,a)@e.T,np.diag([-1,1,1,1]),atol=2e-10,rtol=0)
    return e

def condition(arrays):
    total=sum(arrays);norm=np.linalg.norm(total)
    size=sum(np.linalg.norm(x) for x in arrays)
    # Exclude small/vanishing total components from componentwise statistics.
    mask=abs(total)>1e-6*norm
    comp=sum(abs(x) for x in arrays)/np.maximum(abs(total),1e-300)
    return dict(sum_term_norms=float(size),total_norm=float(norm),
                norm_condition=float(size/max(norm,1e-300)),
                max_component_condition_nonzero=float(max(comp[mask])) if mask.any() else None,
                excluded_small_components=int((~mask).sum()),
                estimated_relative_norm_error_if_each_term_has_1e_minus12_relative_error=float(1e-12*size/max(norm,1e-300)))

def sensitivity(result):
    import lorenz_metric as metric
    from lorenz_jet import Jet
    original_data=metric._homogeneous_radial_data
    original_dtype=Jet.coefficient_dtype
    rows=[]
    for a in [0.,.8771530275949366]:
        for ell in [1,18]:
            for r in [19.9995,20.0005]:
                e=frame(r,1.1,a)
                base=e@values(metric.nonstatic_metric(r,1.1,20.,a,ell,1,order=6)[1])@e.T
                partition=e@sum(pieces(r,1.1,20.,a,ell,1).values())@e.T
                scale=np.linalg.norm(base)
                trial=[]
                for seed in [11,29,53]:
                    rng=np.random.default_rng(seed)
                    factors={spin:(rng.uniform(-1,1)+1j*rng.uniform(-1,1),
                                   rng.uniform(-1,1)+1j*rng.uniform(-1,1)) for spin in [-2,-1,0,1,2]}
                    epsilon=1e-12
                    def perturbed(spin,ee,mm,aa,ww,rr,bc):
                        lam,v,dv=original_data(spin,ee,mm,aa,ww,rr,bc)
                        fv,fd=factors[spin]
                        return lam,v*(1+epsilon*fv),dv*(1+epsilon*fd)
                    metric._homogeneous_radial_data=perturbed
                    try:
                        changed=e@values(metric.nonstatic_metric(r,1.1,20.,a,ell,1,order=6)[1])@e.T
                    finally:
                        metric._homogeneous_radial_data=original_data
                    trial.append(dict(seed=seed,epsilon=epsilon,
                          relative_metric_norm_change=float(np.linalg.norm(changed-base)/scale),
                          amplification=float(np.linalg.norm(changed-base)/scale/epsilon)))
                Jet.coefficient_dtype=np.clongdouble
                try:
                    extended=e@values(metric.nonstatic_metric(r,1.1,20.,a,ell,1,order=6)[1])@e.T
                finally:
                    Jet.coefficient_dtype=original_dtype
                row=dict(a=a,ell=ell,r=r,trial=trial,
                         partition_vs_production_relative_norm_change=float(np.linalg.norm(partition-base)/scale),
                         extended_Jet_relative_metric_norm_change=float(np.linalg.norm(extended-base)/scale))
                rows.append(row);print(json.dumps(row),flush=True)
    result['radial_input_sensitivity']=rows
    result['radial_input_sensitivity_limit']='Perturbs value and first derivative independently by <=sqrt(2)*1e-12 inside homogeneous_field_jet, holding sourced amplitudes, trace/kappa radial data, and angular inputs fixed. This probes conditioning, not the actual error of a valid radial solve.'
    result['extended_Jet_limit']='Changes only Taylor coefficient arithmetic to complex long double; pybhpt radial/angular inputs and mass-variation data remain double.'

def main():
    result=dict(status='running_local_reconstruction_conditioning_audit',rows=[],
      limitations=['Finite local sample, individual ell=1,2,18; not a summed ell<=18 metric or final flux.',
        'Term norm condition numbers do not bound internal cancellation inside each Taylor differential operator.',
        'Componentwise ratios exclude components smaller than 1e-6 of total Frobenius norm.',
        'The error estimate assumes independent 1e-12 relative errors of already reconstructed terms, not unprocessed radial data.'],
      source_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()
        for f in [Path(__file__),ROOT/'src/lorenz_metric.py',ROOT/'src/lorenz_jet.py',ROOT/'src/lorenz_kappa.py']})
    def save():
        tmp=OUT.with_suffix('.tmp');tmp.write_text(json.dumps(result,indent=2)+'\n');tmp.replace(OUT)
    for a in [0.,.8771530275949366]:
        rp=1+np.sqrt(1-a*a)
        for ell in [1,2,18]:
            for r in [rp+.001,3.,19.9995,20.0005,40.]:
                t=time.perf_counter();parts=pieces(r,1.1,20.,a,ell,1);e=frame(r,1.1,a)
                physical={k:e@h@e.T for k,h in parts.items()}
                sectors={prefix:sum(h for k,h in physical.items() if k.startswith(prefix))
                         for prefix in ['s0','s1']+(['s2'] if ell>=2 else [])}
                row=dict(a=a,ell=ell,m=1,r=r,theta=1.1,order=6,
                   term_norms={k:float(np.linalg.norm(h)) for k,h in physical.items()},
                   all_terms=condition(list(physical.values())),
                   sector_sum=condition(list(sectors.values())),
                   scalar_parts=condition([h for k,h in physical.items() if k.startswith('s0')]),
                   tensor_parts=condition([h for k,h in physical.items() if k.startswith('s2')]) if ell>=2 else None,
                   seconds=time.perf_counter()-t)
                result['rows'].append(row);save()
                print(json.dumps(dict(a=a,ell=ell,r=r,condition=row['all_terms']['norm_condition'],
                    component=row['all_terms']['max_component_condition_nonzero'],seconds=row['seconds'])),flush=True)
    sensitivity(result)
    result['status']='completed_local_reconstruction_conditioning_audit';save()

if __name__=='__main__':main()
