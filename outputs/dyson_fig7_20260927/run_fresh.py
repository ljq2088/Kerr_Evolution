"""Independent resumable Dyson Fig.7 finite-resolution diagnostic (88 channels)."""
import os,sys,json,hashlib,traceback,argparse
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from report_environment_forced_mode import argument_parser,build_cloud,build_metric,source_grid,run
from environment_metric_sampling import precompute_metric
from environment_lorenz_mode import ConjugateMetricMode
from environment_source import angular_mode
from source_provenance import source_fingerprint,validate_saved_samples,local_dependency_hashes
from plot_comparison import plot

def save(path,data):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(data,indent=2,default=str));tmp.replace(path)

def main():
    os.nice(10)
    status={'status':'running','pid':os.getpid(),'started':datetime.now(timezone.utc).isoformat(),'expected_channels':88,'completed':[],'failures':[], 'limitations':['Finite-resolution diagnostic; not a convergence-certified reproduction.','Static metric uses the existing L18 matching candidate, with its SHA recorded.','Metric spheroidal truncation L18 is not claimed equivalent to the paper spherical truncation.']}
    save(OUT/'execution.json',status)
    code_hashes=local_dependency_hashes(ROOT/'src',['report_environment_forced_mode','environment_metric_sampling'])
    matching=ROOT/'docs/environment_reproduction/static_tetrad_a0.877153_r20_L18_q32_j18_free18_paper_eps5e-06.json'
    status['matching_sha256']=hashlib.sha256(matching.read_bytes()).hexdigest()
    status['code_hashes']=code_hashes
    # Process both signs together; share only the metric, never conjugate the scalar cloud.
    cloud=None;rows={}
    for mg in list(range(1,14))+[0]:
        if code_hashes!=local_dependency_hashes(ROOT/'src',['report_environment_forced_mode','environment_metric_sampling']):
            raise RuntimeError('Numerical source code changed during run; stop to avoid mixed provenance')
        plan=[(ell,m) for ell in range(2,13) for m in range(-ell,ell+1,2) if abs(m-1)==mg]
        if not plan:continue
        args=argument_parser().parse_args(['--metric-m',str(mg),'--metric-ellmax','18','--angular-order','40','--radial-order','12','--horizon-order','64','--horizon-log','--source-inner-offset','0.0002','--source-outer-radius','320','--green-horizon-offset','0.0001','--green-outer-radius','4000','--infinity-method','coulomb']+(['--static-matching',str(matching)] if mg==0 else []))
        status['active_metric_m']=mg;save(OUT/'execution.json',status)
        try:
            if cloud is None:cloud=build_cloud(args)
            status['cloud']={'alpha':cloud.mu,'a':cloud.a,'omega':cloud.omega}
            base=build_metric(args,cloud)
            _,radii,_=source_grid(args,cloud)
            theta=np.arccos(np.polynomial.legendre.leggauss(args.angular_order)[0])
            base,audit=precompute_metric(base,radii,theta,ROOT/'outputs/metric_cache',workers=2)
            status.setdefault('metric_audits',{})[str(mg)]=audit
            for ell,m in plan:
                selected=base if m-1==mg else ConjugateMetricMode(base)
                a=argparse.Namespace(**vars(args));a.scalar_ell=ell;a.metric_m=m-1;a.output=OUT/f'mode_l{ell}_m{m}.json'
                path,result=run(a,cloud=cloud,metric=selected)
                validate_saved_samples(result,source_fingerprint())
                p=result['parameters'];points=[v for v in result['radial_response_at_panel_boundaries'] if v['r']==20.]
                if len(points)!=1:raise ValueError('Exact particle panel absent')
                z=complex(*points[0]['field']);ang=angular_mode(np.pi/2,ell,m,cloud.a**2*(p['omega']**2-cloud.mu**2))[0]
                rows[ell,m]={'radial':z,'physical':z*ang,'file':str(path),'sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest()}
                status['completed'].append({'ell':ell,'m':m,'file':str(path)})
                multipoles=[]
                for l in range(2,13):
                    missing=[mm for mm in range(-l,l+1,2) if (l,mm) not in rows]
                    row={'ell':l,'complete':not missing,'missing_m':missing}
                    if not missing:
                        row.update(abs_particle_sum_per_epsilon_q=float(abs(sum(rows[l,mm]['physical'] for mm in range(-l,l+1,2)))/cloud.mu**3),abs_radial_sum_per_epsilon_q=float(abs(sum(rows[l,mm]['radial'] for mm in range(-l,l+1,2)))/cloud.mu**3))
                    multipoles.append(row)
                save(OUT/'fresh_particle_field.json',{'status':'finite_resolution_not_converged','multipoles':multipoles,'inputs':[{k:v for k,v in row.items() if k in ('file','sha256')} for row in rows.values()]})
                plot(OUT/'fresh_particle_field.json','FRESH FINITE-RESOLUTION CALCULATION — convergence pending','fresh_vs_paper')
                save(OUT/'execution.json',status)
        except Exception:
            status['failures'].append({'metric_m':mg,'traceback':traceback.format_exc()});save(OUT/'execution.json',status)
    status['status']='completed_finite_resolution_not_converged' if len(status['completed'])==88 else 'incomplete_with_failures'
    status['finished']=datetime.now(timezone.utc).isoformat();save(OUT/'execution.json',status)
if __name__=='__main__':
    try:main()
    except Exception:
        path=OUT/'execution.json';state=json.loads(path.read_text()) if path.exists() else {};state.update(status='failed',error=traceback.format_exc());save(path,state);raise
