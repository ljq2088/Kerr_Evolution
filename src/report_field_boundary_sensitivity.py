"""Fixed-source threshold response sensitivity to outgoing boundary construction.
This diagnostic intentionally includes unconverged boundary controls; they are
never selected by closeness to the paper. It does not reproduce author code.
"""
import hashlib,json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from environment_radial import RadialGreen
from environment_response import SampledResponse
from environment_source import angular_mode
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'docs/environment_reproduction'
def enc(z):
 z=np.asarray(z,complex);return np.stack([z.real,z.imag],axis=-1).tolist()
def run(job):
 orbit,method,outer=job
 p=D/f'forced_mode_nr8_nt18_L18_alpha0.3_rp{orbit}_mg1_sl2_gh0.0001_go32000_coulomb_inner0.0005_outer320_log_h32.json'
 data=json.loads(p.read_text());v=data['parameters'];a=v['metric']['a'];om=v['omega']
 g=RadialGreen(a,v['alpha'],om,2,2,rmax=outer,offset=v['green_horizon_offset'],rtol=1e-11,infinity_method=method)
 rr=np.array([s['r'] for s in data['samples']]);J=np.array([complex(*s['source']) for s in data['samples']]);ww=np.array([s['weight'] for s in data['samples']])
 response=SampledResponse(g,v['source_panels'],rr,J,log_first=True)
 r=np.array([2.,5.,10.,20.,40.,50.,75.,100.,150.,200.,300.]);R=response.evaluate(r)[0]
 S=angular_mode(np.pi/2,2,2,a*a*(om*om-v['alpha']**2))[0]
 return dict(orbit=orbit,method=method,outer_radius=outer,infinity_series_last_term_relative=g.series_last_term_relative,r=r.tolist(),equatorial_scalar22_per_q_epsilon=enc(R*S/v['alpha']**3),ZH=enc(np.dot(ww,g.upsol.sol(rr)[0]*J)/g.w0),ZI=enc(np.dot(ww,g.insol.sol(rr)[0]*J)/g.w0),wronskian_relative_spread=float(np.max(abs(g.wronskian(np.geomspace(g.rmin,320,80))/g.w0-1))),input_file=str(p.relative_to(ROOT)),input_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def main():
 configs=[('coulomb',32000.),('coulomb',64000.),('coulomb',1000.),('series',1000.),('series',4000.)]
 with ProcessPoolExecutor(4) as pool:rows=list(pool.map(run,[(r,m,b) for r in [41.6,41.8] for m,b in configs]))
 for row in rows:
  base=next(b for b in rows if b['orbit']==row['orbit'] and b['method']=='coulomb' and b['outer_radius']==32000)
  y=np.array(row['equatorial_scalar22_per_q_epsilon']);yb=np.array(base['equatorial_scalar22_per_q_epsilon']);z=y[:,0]+1j*y[:,1];zb=yb[:,0]+1j*yb[:,1]
  row['relative_complex_difference']= (abs(z-zb)/abs(zb)).tolist();row['amplitude_ratio']=(abs(z)/abs(zb)).tolist()
 out=ROOT/'docs/field_alignment_20260917/boundary_sensitivity.json'
 result=dict(status='fixed_historical_source_boundary_control_not_author_error_identification',rows=rows,source_scaled=False,production_solver_changed=False,limitations=['Single scalar22 channel, frozen source; full source and angular functions were not regenerated.','No author outer radius or boundary algorithm is public; matching any failed boundary control would not identify their implementation.','Small Wronskian variation tests integration consistency, not correctness of the imposed outer boundary.'],implementation_sha256={str(q.relative_to(ROOT)):hashlib.sha256(q.read_bytes()).hexdigest() for q in [Path(__file__),ROOT/'src/environment_radial.py',ROOT/'src/environment_response.py']})
 out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 for row in rows:print(row['orbit'],row['method'],row['outer_radius'],row['amplitude_ratio'][7],row['infinity_series_last_term_relative'],flush=True)
if __name__=='__main__':main()
