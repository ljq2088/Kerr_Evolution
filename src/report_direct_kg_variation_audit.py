"""Audit the Lorenz environmental source against finite variation of full KG.

Directly evaluates Box in coordinate divergence form for g+epsilon*h. No
Lorenz condition, Christoffel symbols, or covariant Hessian enters that path.
The same reconstructed h and normalized cloud are physical inputs; this is
an operator/constraint cross-check, not independent reconstruction validation.
"""
import json,hashlib
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np


def enc(v):
    v=np.asarray(v);return np.stack((v.real,v.imag),axis=-1).tolist()


def kg(metric,derivative,phi,grad,second,mu):
    inv=np.linalg.inv(metric)
    contracted=np.zeros(4,complex)
    for a in range(4):
        dinv=-inv@derivative[a]@inv
        contracted+=dinv[a]+.5*inv[a]*np.trace(inv@derivative[a])
    return np.einsum('ab,ab->',inv,second)+contracted@grad-mu**2*phi


def evaluate(case):
    from environment_source import ThresholdCloud,angular_mode,kerr_metric
    from lorenz_metric import nonstatic_metric
    cloud=ThresholdCloud(alpha=.3)
    r,theta,L,mg=case;r=cloud.rp+.005 if r=='near_horizon' else float(r)
    h=np.zeros((4,4),complex);dh=np.zeros((4,4,4),complex)
    for ell in range(1,L+1):
        g,piece=nonstatic_metric(r,theta,20.,cloud.a,ell,1,order=8)
        h+=np.array([[v.value for v in row] for row in piece])
        dh+=np.array([[[g.partial(v,b).value for v in row] for row in piece] for b in range(4)])
    if mg==-1:h,dh=h.conjugate(),dh.conjugate()
    g0=kerr_metric(r,theta,cloud.a);dg=np.zeros((4,4,4),complex)
    dg[1]=kerr_metric(r+1e-25j,theta,cloud.a).imag/1e-25
    dg[2]=kerr_metric(r,theta+1e-25j,cloud.a).imag/1e-25
    R,Rp=cloud.radial(r)[:,0];Rpp=cloud.rhs(r,[R,Rp])[1]
    S,Sp,A=angular_mode(theta,1,1,cloud.c2)
    Spp=-np.cos(theta)/np.sin(theta)*Sp-(A+cloud.c2*np.cos(theta)**2-1/np.sin(theta)**2)*S
    phi=R*S;grad=np.array([-1j*cloud.omega*phi,Rp*S,R*Sp,1j*phi])
    second=np.empty((4,4),complex)
    second[0,:]=second[:,0]=-1j*cloud.omega*grad
    second[3,:]=second[:,3]=1j*grad
    second[1,1]=Rpp*S;second[2,2]=R*Spp;second[1,2]=second[2,1]=Rp*Sp
    reference=cloud.lorenz_source(r,theta,h)
    scale=max(float(np.linalg.norm(np.linalg.inv(g0)@h,ord=2)),1e-30)
    rows=[]
    for q in (.01,.005,.0025,.00125,.000625,.0003125):
        eps=q/scale
        estimate=-(kg(g0+eps*h,dg+eps*dh,phi,grad,second,cloud.mu)-kg(g0-eps*h,dg-eps*dh,phi,grad,second,cloud.mu))/(2*eps)
        rows.append(dict(relative_metric_step=q,epsilon=eps,source=enc(estimate),relative_difference=float(abs(estimate-reference)/max(abs(reference),1e-300))))
    last=complex(*rows[-1]['source']);prev=complex(*rows[-2]['source']);extrap=(4*last-prev)/3
    prevex=(4*prev-complex(*rows[-3]['source']))/3
    result=dict(r=r,theta=theta,metric_ellmax=L,m_g=mg,jet_order=8,source_lorenz=enc(reference),source_direct_extrapolated=enc(extrap),relative_difference=float(abs(extrap-reference)/max(abs(reference),1e-300)),relative_extrapolation_change=float(abs(extrap-prevex)/max(abs(extrap),1e-300)),background_KG_residual=enc(kg(g0,dg,phi,grad,second,cloud.mu)),metric_perturbation_scale=scale,rows=rows)
    print('completed',case,'relative difference',result['relative_difference'],flush=True)
    return result


def main():
    root=Path(__file__).resolve().parents[1];out=root/'docs/environment_reproduction/direct_kg_variation_audit_20260916.json'
    cases=[(2.,1.1,4,-1),(10.,.7,4,-1),(10.,1.2,4,-1),(19.9,1.1,4,-1),(20.1,1.1,4,-1),(30.,.7,4,-1),('near_horizon',1.1,4,-1),(10.,1.2,4,1)]
    names=['report_direct_kg_variation_audit.py','environment_source.py','lorenz_metric.py','lorenz_ghp.py','lorenz_kappa.py','lorenz_spin1_chiral.py']
    result=dict(status='running',parameters=dict(alpha=.3,rp=20.),implementation_sha256={n:hashlib.sha256((root/'src'/n).read_bytes()).hexdigest() for n in names},rows=[],limitations=['Same reconstructed metric and cloud inputs; not an independent physical source or normalization.','Only sampled points and ellmax4, away from the orbital surface; distributional terms require separate checks.','Finite epsilon extrapolation of analytic complex metric perturbation; no fitted scales or phases.'])
    with ProcessPoolExecutor(max_workers=2) as pool:
        jobs=[pool.submit(evaluate,c) for c in cases]
        for job in as_completed(jobs):
            result['rows'].append(job.result());out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    result['status']='completed_direct_full_KG_variation_not_total_flux_validation';result['maximum_relative_difference']=max(x['relative_difference'] for x in result['rows']);out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
