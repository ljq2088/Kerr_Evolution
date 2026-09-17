"""Counterfactual omitted Appendix term: unchanged correct production h:H source.
Read cached Lorenz metric, project only the analytically identified difference,
and propagate with exactly the baseline Green function. This is NOT evidence
that the author's program used the printed expression.
"""
from pathlib import Path
import sys,json,hashlib,time
import numpy as np
from scipy.optimize import minimize_scalar
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from environment_source import ThresholdCloud,angular_mode,project_source
from environment_lorenz_mode import LorenzMetricMode
from environment_metric_sampling import precompute_metric
from environment_dense_metric import configure_metric_backend
from environment_response import SampledResponse
from environment_cloud import mode_flux
from source_provenance import source_fingerprint,validate_saved_samples
BASE=Path(__file__).resolve().parent

def enc(z):return [float(np.real(z)),float(np.imag(z))]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def minima(response,ang):
    r=np.linspace(5.,200.,1951);v=abs(response.evaluate(r)[0]*ang)
    indices=np.where((v[1:-1]<v[:-2])&(v[1:-1]<v[2:]))[0]+1;rows=[]
    for i in indices:
        opt=minimize_scalar(lambda z:float(abs(response.evaluate([z])[0][0]*ang)**2),bounds=(r[i-1],r[i+1]),method='bounded',options={'xatol':1e-9})
        rows.append(dict(r=float(opt.x),amplitude=float(np.sqrt(opt.fun))))
    return rows

def run(orbit,cloud,fp):
    begin=time.perf_counter();tag=str(orbit).replace('.','p')
    path=ROOT/f'docs/field_alignment_20260917/nonstatic/rp{tag}_sl2_sm2_L6_q12_nr8_h32.json'
    old=json.loads(path.read_text());validate_saved_samples(old,fp);p=old['parameters']
    a=cloud.a;mu=cloud.mu;w=p['omega'];radii=np.array([v['r'] for v in old['samples']]);weights=np.array([v['weight'] for v in old['samples']]);J=np.array([complex(*v['source']) for v in old['samples']])
    if abs(a-p['metric']['a'])>1e-14 or p['angular_order']!=12 or p['metric']['ellmax']!=6:raise ValueError('Baseline/cloud mismatch')
    x,wt=np.polynomial.legendre.leggauss(12);theta=np.arccos(x)
    original_metric=LorenzMetricMode(orbit,a,1,6)
    cached,audit=precompute_metric(original_metric,radii,theta,ROOT/'outputs/metric_cache',workers=4)
    print('CACHE',orbit,audit,flush=True)
    _,fresh=project_source(cloud,radii,orbit,2,2,cached,ntheta=12)
    check=float(np.linalg.norm(fresh-J)/np.linalg.norm(J))
    if check>1e-10:raise ValueError(f'Baseline source reconstruction mismatch {check}')
    Sc=angular_mode(theta,1,1,cloud.c2)[0];Sout=angular_mode(theta,2,2,a*a*(w*w-mu*mu))[0]
    lamc=cloud.lam+a*a*cloud.omega**2-2*a*cloud.omega
    missing=[]
    for r in radii:
        delta=r*r-2*r+a*a;sig=r*r+a*a*x*x
        lp=np.array([(r*r+a*a)/delta,1.,0.,a/delta]);lm=np.array([-(r*r+a*a)/delta,1.,0.,-a/delta])
        hlplm=np.array([lp@cached(float(r),float(t))@lm for t in theta])
        Phi=cloud.radial(r)[0,0]*Sc
        difference=delta*hlplm*(mu*mu*r*r+lamc)*Phi/(2*sig)
        missing.append(2*np.pi*np.dot(wt,Sout*difference))
    missing=np.asarray(missing)
    baseline=SampledResponse.from_report(old);green=baseline.green
    args=(green,p['source_panels'],radii)
    delta_response=SampledResponse(*args,missing,log_first=True)
    printed=SampledResponse(*args,J-missing,log_first=True)
    radii_fixed=np.array([2.,3.,5.,10.,20.,30.,40.,50.,70.,100.,150.,180.,200.]);dense=np.linspace(5.,200.,1951)
    ang=angular_mode(np.pi/2,2,2,a*a*(w*w-mu*mu))[0]/mu**3
    b=baseline.evaluate(radii_fixed)[0]*ang;d=delta_response.evaluate(radii_fixed)[0]*ang;t=printed.evaluate(radii_fixed)[0]*ang
    bd=baseline.evaluate(dense)[0]*ang;td=printed.evaluate(dense)[0]*ang
    z={}
    for name,attr in [('up_coefficient','up_coefficient'),('z_h','horizon_coefficient')]:
        zb,zd,zt=[getattr(v,attr) for v in [baseline,delta_response,printed]]
        z[name]=dict(correct_direct=enc(zb),omitted_contribution=enc(zd),counterfactual_printed=enc(zt),relative_complex_change=float(abs(zt-zb)/max(abs(zb),1e-300)),linearity_absolute_error=float(abs(zb-zd-zt)))
    probes=np.array([radii[0],3.,20.,orbit,100.,320.])
    result=dict(status='counterfactual_published_appendix_omission_control_complete',orbital_radius=orbit,parameters=p,
        baseline_path=str(path.relative_to(ROOT)),baseline_sha256=sha(path),source_provenance=fp,
        implementation_sha256=sha(__file__),metric_cache=audit,baseline_source_reconstruction_relative_L2=check,
        formula=dict(target='Difference in projected J=integral S22* Sigma source dOmega',
            local_difference='Delta*h_lplus_lminus*(mu^2*r^2+lambda_cloud)*Phi_cloud/(2*Sigma)',lambda_cloud=float(lamc)),
        source_difference_relative_L2=float(np.linalg.norm(missing)/np.linalg.norm(J)),
        source_difference_weighted_relative_L2=float(np.sqrt(np.sum(weights*abs(missing)**2)/np.sum(weights*abs(J)**2))),
        source_samples=[dict(r=float(r),weight=float(wt),correct_direct=enc(j),omitted=enc(q),counterfactual_printed=enc(j-q)) for r,wt,j,q in zip(radii,weights,J,missing)],
        amplitudes=z,infinity_propagating=bool(green.propagating),
        z_inf=dict(correct_direct=enc(baseline.up_coefficient if green.propagating else 0j),counterfactual_printed=enc(printed.up_coefficient if green.propagating else 0j)),
        fluxes={label:mode_flux(w,2,cloud.omega,1,mu,a,response.up_coefficient,response.horizon_coefficient) for label,response in [('correct_direct',baseline),('counterfactual_printed',printed)]},
        fixed_radius_fields=[dict(r=float(r),correct_direct=enc(zb),omitted_contribution=enc(zd),counterfactual_printed=enc(zt),relative_complex_change=float(abs(zt-zb)/max(abs(zb),1e-300)),amplitude_ratio=float(abs(zt)/max(abs(zb),1e-300))) for r,zb,zd,zt in zip(radii_fixed,b,d,t)],
        dense_profile=dict(radii=dense.tolist(),correct_amplitude=abs(bd).tolist(),counterfactual_amplitude=abs(td).tolist(),
            relative_complex_L2=float(np.linalg.norm(td-bd)/np.linalg.norm(bd)),relative_amplitude_L2=float(np.linalg.norm(abs(td)-abs(bd))/np.linalg.norm(abs(bd))),
            correct_radial_amplitude_minima=minima(baseline,ang),counterfactual_radial_amplitude_minima=minima(printed,ang)),
        linearity_relative_error=float(np.linalg.norm(b-d-t)/np.linalg.norm(b)),
        green_wronskian_relative_spread=float(np.max(abs(green.wronskian(probes)/green.w0-1))),
        elapsed_seconds=time.perf_counter()-begin,
        interpretation='Correct production source is direct h:H and already contains this term. Counterfactual deliberately removes it to measure a printed-formula discrepancy, not to align a figure.',
        limitations=['Finite metric L6/angular12/radial8/horizon32 baseline; no new convergence claim.',
            'Exact synchronized spin differs from Li fixed a=.88; all other source and boundary settings unchanged.',
            'No evidence of which source formula the authors executed; do not attribute the counterfactual to their code.',
            'Only scalar22 response; minima are minima of its radial absolute amplitude, not full complex-field zeros.'])
    if source_fingerprint()!=fp:raise RuntimeError('Production changed during diagnostic')
    out=BASE/f'li_appendix_omission_rp{tag}.json';out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('COMPLETE',orbit,'sourceL2',result['source_difference_relative_L2'],'fieldL2',result['dense_profile']['relative_complex_L2'],'ZHchange',z['z_h']['relative_complex_change'],flush=True)
    return result

def main():
    configure_metric_backend(False);fp=source_fingerprint();cloud=ThresholdCloud(alpha=.3)
    rows=[run(orbit,cloud,fp) for orbit in [42.1,41.1]]
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for axis,row in zip(ax,rows):
        q=row['dense_profile'];axis.plot(q['radii'],q['correct_amplitude'],label='Direct covariant source (correct)');axis.plot(q['radii'],q['counterfactual_amplitude'],ls='--',label='Printed Appendix counterfactual');axis.set(xlabel='r/M',ylabel=r'$|R_{22}S_{22}(\pi/2)|/\alpha^3$',title=f"Orbit {row['orbital_radius']} M",yscale='log');axis.legend(fontsize=8);axis.grid(alpha=.2)
    fig.savefig(BASE/'li_appendix_omission_control.png',dpi=170)
if __name__=='__main__':main()
