"""Warm only mg5/6 metric caches for the independent Li-field orchestrator.

Uses its unchanged argument/grid factory. Never writes scalar response files.
"""
import argparse,json,os,time
from pathlib import Path
import numpy as np
from report_li_field_alignment_runs import make_args,orbit_tag,stamp,save,sha,BASE,ROOT

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--rp',type=float,required=True);parser.add_argument('--workers',type=int,default=3);args=parser.parse_args()
 from environment_source import ThresholdCloud
 from environment_lorenz_mode import LorenzMetricMode
 from environment_metric_sampling import precompute_metric
 from environment_dense_metric import configure_metric_backend
 from report_environment_forced_mode import source_grid
 from source_provenance import source_fingerprint
 configure_metric_backend(False);cloud=ThresholdCloud(alpha=.3);fp=source_fingerprint()
 out=BASE/f'metric_warm_rp{orbit_tag(args.rp)}_manifest.json';main_manifest=BASE/f'nonstatic_rp{orbit_tag(args.rp)}_manifest.json'
 result=dict(status='warming_only_metric_cache',started_utc=stamp(),pid=os.getpid(),orbital_radius=args.rp,workers=args.workers,metric_ms=[5,6],source_provenance=fp,implementation_sha256=sha(__file__),orchestrator_sha256=sha(ROOT/'src/report_li_field_alignment_runs.py'),groups=[],scalar_responses_written=False)
 save(out,result)
 for mg,ell,m in [(5,4,-4),(6,5,-5)]:
  state=json.loads(main_manifest.read_text())
  if state.get('active_metric_m',0)>=mg and not any(v['metric_m']==mg for v in state.get('groups',[])):
   result.update(status='stopped_main_sampler_already_owns_target_group',conflicting_metric_m=mg);save(out,result);raise RuntimeError(f'Main sampler already active at mg={state.get("active_metric_m")}; do not co-write target cache')
  settings=make_args(args.rp,ell,m);_,radii,_=source_grid(settings,cloud);theta=np.arccos(np.polynomial.legendre.leggauss(12)[0]);base=LorenzMetricMode(args.rp,cloud.a,mg,6)
  result.update(active_metric_m=mg,active_metric_started_utc=stamp());save(out,result)
  print(f'WARM ORBIT {args.rp} MG {mg} START {len(radii)} RADII',flush=True)
  sampled,audit=precompute_metric(base,radii,theta,ROOT/'outputs/metric_cache',args.workers)
  after=source_fingerprint()
  if after!=fp:raise RuntimeError('Production source changed during metric cache warming')
  result['groups'].append(dict(metric_m=mg,completed_utc=stamp(),precomputation=audit,source_provenance_verified=True,grid_args=dict(ell=ell,m=m,source_outer_radius=settings.source_outer_radius,scalar_output_path_untouched=str(settings.output)),radii=radii.tolist(),theta=theta.tolist(),metric_provenance=base.provenance))
  save(out,result);print(f'WARM ORBIT {args.rp} MG {mg} COMPLETE '+json.dumps(audit),flush=True)
 result.update(status='metric_cache_warming_complete',completed_utc=stamp(),source_provenance_after=source_fingerprint());result.pop('active_metric_m',None);save(out,result)
 print(out,flush=True)
if __name__=='__main__':main()
