"""Propagate saved author/local metric projections through one shared scalar source.

This is NOT original-author environment-source code. Both inputs are truncated
identically to protected output spherical j=1..3; neither J is the full source.
No metric is regenerated. A real metric's m_g=-1 coefficient is the complex
conjugate of its m_g=+1 BL tensor at phi=t=0, without a factor of two.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from pybhpt.swsh import Yslm
from paper_full_tetrad import SPINS,to_bl,project_metric
from environment_source import ThresholdCloud,angular_mode
ROOT=Path(__file__).resolve().parents[1]

def enc(z):
    z=np.asarray(z,complex);return np.stack([z.real,z.imag],axis=-1).tolist()
def dec(z):
    z=np.asarray(z);return z[...,0]+1j*z[...,1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def unweight(q,r,t,a):
    v=np.array(q,complex);rho=r+1j*a*np.cos(t);sigma=abs(rho)**2;delta=r*r-2*r+a*a
    v[4]/=rho;v[5]/=rho.conjugate();v[6]/=rho.conjugate();v[7]/=rho
    v[8]/=sigma*delta
    v[9]=sigma*q[9]-delta*v[8]
    return v

def project_one(cloud,r,coeff,m,jcut,nquad):
    x,w=np.polynomial.legendre.leggauss(nquad);theta=np.arccos(x)
    omega=cloud.omega-m/(20**1.5+cloud.a)
    angular_out=angular_mode(theta,0,0,cloud.a**2*(omega**2-cloud.mu**2))[0]
    J_by_component=np.zeros(10,complex);source_integrand=[];metric_norm=[];termabs_integral=0.;roundtrip=0.
    for k,t in enumerate(theta):
        q=np.array([sum(coeff[j,c]*Yslm(s,j,m,t) for j in range(abs(m),jcut+1) if j>=abs(s)) for c,s in enumerate(SPINS)])
        _,hessian,inverse=cloud.hessian(r,t)
        H=to_bl(unweight(q,r,t,cloud.a),r,t,cloud.a)
        roundtrip=max(roundtrip,float(np.max(abs(project_metric(H,r,t,cloud.a)-unweight(q,r,t,cloud.a)))))
        hcontra=inverse@H.conjugate()@inverse
        contraction_terms=hcontra*hessian
        factor=2*np.pi*w[k]*angular_out[k]*(r*r+cloud.a**2*x[k]*x[k])
        source_integrand.append(np.sum(contraction_terms));metric_norm.append(float(np.linalg.norm(H)))
        termabs_integral+=abs(factor)*np.sum(abs(contraction_terms))
        for c in range(10):
            v=np.zeros(10,complex);v[c]=q[c]
            h=to_bl(unweight(v,r,t,cloud.a),r,t,cloud.a).conjugate()
            J_by_component[c]+=factor*np.sum((inverse@h@inverse)*hessian)
    integrand=np.asarray(source_integrand)
    factors=2*np.pi*w*angular_out*(r*r+cloud.a**2*x*x)
    J=np.dot(factors,integrand)
    scale=abs(J)
    return dict(J=enc(J),J_complex=J,J_by_weighted_component=enc(J_by_component),component_sum_consistency_abs=float(abs(sum(J_by_component)-J)),weighted_component_cancellation=float(np.sum(abs(J_by_component))/scale) if scale else None,angular_cancellation=float(np.dot(abs(factors),abs(integrand))/scale) if scale else None,coordinate_contraction_and_angular_cancellation=float(termabs_integral/scale) if scale else None,tetrad_roundtrip_max_abs=roundtrip,integrand_complex=integrand,metric_norm=metric_norm)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--comparison',required=True);parser.add_argument('--output',required=True);parser.add_argument('--quadratures',type=int,nargs='+',default=[48,96]);args=parser.parse_args()
    path=Path(args.comparison);original=json.loads(path.read_text());protected=original['protected_output_ells']
    if protected!=[1,2,3] or original['m']!=1 or original['rp']!=20:raise ValueError('This bounded diagnostic requires actual m1/rp20/protected j1..3 inputs')
    start=time.perf_counter();cloud=ThresholdCloud(alpha=.3)
    if abs(cloud.a-original['a'])>1e-12:raise ValueError('Cloud spin does not match author metric')
    result=dict(status='completed_author_metric_input_propagation',classification='original_author_metric_input_shared_local_environment_source',scope='No original author cloud/Hessian/environment-source implementation is run. Both author and local saved metrics truncated to output spherical j1..3, same cloud/source/quadrature. Neither J is the full environmental source or flux.',input_comparison=str(path),input_comparison_sha256=sha(path),raw_author_path=original['reference_path'],raw_author_sha256=original['reference_sha256'],reference_json_normalization=original.get('reference_json_normalization'),parameters=dict(a=cloud.a,alpha=cloud.mu,omega_cloud=cloud.omega,cloud_mass=cloud.mass,rp=20,metric_m_input=1,metric_m_after_conjugation=-1,scalar_ell=0,scalar_m=0,protected_output_j=protected,quadratures=args.quadratures),formula='J_00=2*pi integral S_00(theta)*Sigma*conj(h_+1)^ab*nabla_a nabla_b Phi_211 sin(theta)dtheta; T23=Sigma*q9-Delta*T01, T01=q8/(Sigma*Delta) using zero-based q',implementation_sha256={n:sha(ROOT/'src'/n) for n in ['report_author_metric_source_propagation.py','paper_full_tetrad.py','environment_source.py','environment_cloud.py']},rows=[])
    for row in original['rows']:
        radius=row['r'];ref=dec(row['reference']);local=dec(row['local'])
        for q in args.quadratures:
            a=project_one(cloud,radius,ref,1,3,q);b=project_one(cloud,radius,local,1,3,q)
            Ja=a.pop('J_complex');Jb=b.pop('J_complex');sa=a.pop('integrand_complex');sb=b.pop('integrand_complex')
            metricdiff=np.linalg.norm((local-ref)[1:4]);metricscale=np.linalg.norm(ref[1:4]);relJ=abs(Jb-Ja)/abs(Ja) if abs(Ja)>0 else None
            out=dict(r=radius,side=row['side'],ntheta=q,author=a,local=b,difference=enc(Jb-Ja),relative_J_difference=float(relJ) if relJ is not None else None,coefficient_euclidean_relative_difference=float(metricdiff/metricscale),source_integrand_l2_relative_difference=float(np.linalg.norm(sb-sa)/np.linalg.norm(sa)),warning='Coefficient Euclidean norm combines differently weighted components; diagnostic only, not an invariant metric norm.')
            result['rows'].append(out)
            print('r',radius,'q',q,'J author',Ja,'local',Jb,'rel',relJ,flush=True)
    result['elapsed_seconds']=time.perf_counter()-start
    result['quadrature_checks']=[]
    for r in original['rows']:
        entries=[e for e in result['rows'] if e['r']==r['r']]
        if len(entries)>1:
            first,last=entries[0],entries[-1]
            result['quadrature_checks'].append(dict(r=r['r'],author_relative=float(abs(dec(first['author']['J'])-dec(last['author']['J']))/abs(dec(last['author']['J']))),local_relative=float(abs(dec(first['local']['J'])-dec(last['local']['J']))/abs(dec(last['local']['J'])))))
    Path(args.output).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
