"""Fresh source-resolution control of Li scalar22; no production modifications."""
import argparse,json,os,time,traceback
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from numpy.polynomial.legendre import legfit,legval
from report_li_field_alignment_runs import make_args,orbit_tag,stamp,save,sha,ROOT
BASE=ROOT/'docs/root_cause_followup_20260917'
def enc(z):
 z=complex(z);return [float(z.real),float(z.imag)]
def relative(new,old):return float(abs(new-old)/max(abs(old),1e-300))
def samples(data):
 rows=data['samples'];return np.array([x['r'] for x in rows]),np.array([complex(*x['source']) for x in rows]),np.array([x['weight'] for x in rows])
def compare(old,new,oldpath,newpath,out):
 from environment_response import SampledResponse
 from environment_source import angular_mode
 from source_provenance import source_fingerprint,validate_saved_samples
 fp=source_fingerprint();validate_saved_samples(old,fp);validate_saved_samples(new,fp)
 p=old['parameters'];q=new['parameters'];source_error=[];rl,jl,wl=samples(old);rh,jh,wh=samples(new)
 hp=1+np.sqrt(1-p['metric']['a']**2);panels=p['source_panels']
 for i,(lo,hi) in enumerate(zip(panels[:-1],panels[1:])):
  ml=(rl>lo)&(rl<hi);mh=(rh>lo)&(rh<hi)
  trans=(lambda x:np.log(x-hp)) if i==0 and p.get('horizon_log_first_panel') else (lambda x:x)
  a,b=trans(np.array([lo,hi]));x=2*(trans(rl[ml])-a)/(b-a)-1;y=2*(trans(rh[mh])-a)/(b-a)-1
  pred=legval(y,legfit(x,jl[ml],len(x)-1));err=pred-jh[mh]
  source_error.append(dict(panel=[lo,hi],coarse_count=int(ml.sum()),fine_count=int(mh.sum()),relative_weighted_L2=float(np.sqrt(np.sum(wh[mh]*abs(err)**2)/np.sum(wh[mh]*abs(jh[mh])**2))),absolute_max=float(max(abs(err))),source_absolute_max=float(max(abs(jh[mh])))))
 low=SampledResponse.from_report(old);high=SampledResponse.from_report(new)
 radii=np.array([2.,3.,5.,10.,20.,30.,40.,50.,70.,100.,150.,180.,200.])
 sl=angular_mode(np.pi/2,2,2,p['metric']['a']**2*(p['omega']**2-p['alpha']**2))[0]/p['alpha']**3
 sh=angular_mode(np.pi/2,2,2,q['metric']['a']**2*(q['omega']**2-q['alpha']**2))[0]/q['alpha']**3
 vl=low.evaluate(radii)[0]*sl;vh=high.evaluate(radii)[0]*sh
 common=[dict(r=float(r),coarse_field=enc(l),fine_field=enc(h),relative_complex_difference=relative(h,l),coarse_amplitude=float(abs(l)),fine_amplitude=float(abs(h))) for r,l,h in zip(radii,vl,vh)]
 dense=np.linspace(5,200,1951);dl=low.evaluate(dense)[0]*sl;dh=high.evaluate(dense)[0]*sh
 def minima(response,ang):
  values=abs(response.evaluate(dense)[0]*ang);indices=np.flatnonzero((values[1:-1]<values[:-2])&(values[1:-1]<values[2:]))+1;result=[]
  for i in indices:
   opt=minimize_scalar(lambda r:float(abs(response.evaluate([r])[0][0]*ang)**2),bounds=(dense[i-1],dense[i+1]),method='bounded',options={'xatol':1e-9})
   result.append(dict(r=float(opt.x),amplitude=float(np.sqrt(opt.fun))))
  return result
 amplitudes={}
 for name in ['z_h','z_inf']:
  a,b=complex(*old[name]),complex(*new[name]);amplitudes[name]=dict(coarse=enc(a),fine=enc(b),relative_complex_difference=relative(b,a))
 for name,attr in [('cumulative_z_h','horizon_coefficient'),('cumulative_z_inf','up_coefficient')]:
  a,b=getattr(low,attr),getattr(high,attr);amplitudes[name]=dict(coarse=enc(a),fine=enc(b),relative_complex_difference=relative(b,a))
 result=dict(status='fresh_source_resolution_control_complete',created_utc=stamp(),source_provenance=fp,inputs_sha256={str(oldpath.relative_to(ROOT)):sha(oldpath),str(newpath.relative_to(ROOT)):sha(newpath)},implementation_sha256=sha(__file__),parameters=dict(coarse=p,fine=q),fixed_radius_complex_fields=common,amplitudes=amplitudes,flux=dict(coarse=old['flux'],fine=new['flux']),source_interpolation_test=dict(definition='Evaluate the original piecewise source polynomial at newly recomputed fine Gauss source nodes; no amplitude or phase fit.',panels=source_error),dense_profile=dict(r_interval=[5,200],points=len(dense),relative_complex_L2=float(np.linalg.norm(dh-dl)/np.linalg.norm(dl)),relative_amplitude_L2=float(np.linalg.norm(abs(dh)-abs(dl))/np.linalg.norm(abs(dl))),coarse_radial_amplitude_minima=minima(low,sl),fine_radial_amplitude_minima=minima(high,sh)),limitations=['Single scalar22 response at fixed synchronized background; not full-mode field or exact a=.88.', 'Radial amplitude minima are minima of |R22|, not exact zeros of the full coherent field.', 'Finer source sampling tests radial discretization only when metric L and angular order are unchanged.', 'Source support and homogeneous boundary prescription retained; this is not a source-cutoff or Green-boundary convergence test.'])
 save(out,result);print('COMPARE_COMPLETE',out,'FIELD_L2',result['dense_profile']['relative_complex_L2'],flush=True)
def main():
 pa=argparse.ArgumentParser();pa.add_argument('--rp',type=float,required=True);pa.add_argument('--workers',type=int,default=4);pa.add_argument('--nr',type=int,default=16);pa.add_argument('--h',type=int,default=64);pa.add_argument('--L',type=int,default=6);pa.add_argument('--q',type=int,default=12);pa.add_argument('--compare-only',action='store_true');a=pa.parse_args()
 from environment_source import ThresholdCloud
 from environment_lorenz_mode import LorenzMetricMode
 from environment_metric_sampling import precompute_metric
 from environment_dense_metric import configure_metric_backend
 from report_environment_forced_mode import source_grid,run
 from source_provenance import source_fingerprint,validate_saved_samples
 configure_metric_backend(False);cloud=ThresholdCloud(alpha=.3);fp=source_fingerprint();BASE.mkdir(parents=True,exist_ok=True)
 tag=f'rp{orbit_tag(a.rp)}_sl2_sm2_L{a.L}_q{a.q}_nr{a.nr}_h{a.h}'
 out=BASE/(tag+'.json');manifest=BASE/(tag+'_manifest.json');comparison=BASE/(tag+'_comparison.json')
 oldpath=ROOT/'docs/field_alignment_20260917/nonstatic'/f'rp{orbit_tag(a.rp)}_sl2_sm2_L6_q12_nr8_h32.json';old=json.loads(oldpath.read_text());validate_saved_samples(old,fp)
 record=dict(status='starting',started_utc=stamp(),pid=os.getpid(),implementation_sha256=sha(__file__),source_provenance=fp,baseline=str(oldpath.relative_to(ROOT)),baseline_sha256=sha(oldpath),parameters=vars(a),output=str(out.relative_to(ROOT)));save(manifest,record)
 started=time.perf_counter()
 try:
  if a.compare_only:
   new=json.loads(out.read_text());validate_saved_samples(new,fp)
  else:
   args=make_args(a.rp,2,2);args.output=out;args.radial_order=a.nr;args.horizon_order=a.h;args.angular_order=a.q;args.metric_ellmax=a.L
   panels,radii,weights=source_grid(args,cloud);theta=np.arccos(np.polynomial.legendre.leggauss(a.q)[0]);metric=LorenzMetricMode(a.rp,cloud.a,1,a.L)
   record.update(status='metric_precompute',source_nodes=len(radii),radii=radii.tolist(),theta=theta.tolist());save(manifest,record)
   print('FRESH_METRIC_START',tag,len(radii),flush=True)
   cached,audit=precompute_metric(metric,radii,theta,ROOT/'outputs/metric_cache',a.workers)
   if source_fingerprint()!=fp:raise RuntimeError('Production source changed during sampling')
   record.update(status='source_projection',metric_sampling=audit);save(manifest,record)
   print('FRESH_SOURCE_START',tag,flush=True);_,new=run(args,cloud=cloud,metric=cached);validate_saved_samples(new,fp)
  record.update(status='response_comparison',source_sha256=sha(out));save(manifest,record)
  compare(old,new,oldpath,out,comparison)
  if source_fingerprint()!=fp:raise RuntimeError('Production source changed during comparison')
  record.update(status='complete',completed_utc=stamp(),elapsed_seconds=time.perf_counter()-started,comparison=str(comparison.relative_to(ROOT)),comparison_sha256=sha(comparison));save(manifest,record)
 except Exception as exc:
  record.update(status='failed_resumable',exception=repr(exc),traceback=traceback.format_exc(),elapsed_seconds=time.perf_counter()-started);save(manifest,record);raise
if __name__=='__main__':main()
