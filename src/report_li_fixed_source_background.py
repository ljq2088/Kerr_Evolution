"""Change the radial propagator/background frequency while holding sources fixed.

This is a sensitivity experiment, never a consistent new-background solution.
The original logarithmic source coordinate is kept for interpolating its samples.
"""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from environment_radial import RadialGreen
from environment_response import SampledResponse
from environment_source import angular_mode
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'docs/field_alignment_20260917';D=ROOT/'docs/environment_reproduction'
POS={(2,2),(3,3),(4,2),(4,4),(5,3),(5,5)}
class OriginalSourceCoordinate:
 def __init__(self,physical,source_rp):self.physical=physical;self.rp=source_rp
 def __getattr__(self,name):return getattr(self.physical,name)
def main():
 pa=argparse.ArgumentParser();pa.add_argument('--omega-cloud',type=float,required=True);pa.add_argument('--spin',type=float,default=.88);pa.add_argument('--output',type=Path,required=True);args=pa.parse_args()
 ref=np.load(D/'li_field_polar_reference_20260917.npz');r=ref['radii'];phi=ref['phi'];rows=[];inputs={};mode_details=[]
 for io,orbit in enumerate([41.1,42.1]):
  values={'all':np.zeros((len(r),len(phi)),complex),'positive':np.zeros((len(r),len(phi)),complex)}
  files=sorted((OUT/'nonstatic').glob(f'rp{str(orbit).replace(".","p")}_sl*_L6_q12_nr8_h32.json'))
  files += [D/f'forced_mode_nr8_nt12_L6_alpha0.3_rp{orbit}_mg0_sl{ell}_gh0.0001_go4000_coulomb_inner0.0005_outer320_log_h32.json' for ell in [3,5]]
  if len(files)!=18:raise ValueError('Require all 18 fixed-source channels')
  for path in files:
   data=json.loads(path.read_text());inputs[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest();p=data['parameters'];l,m=p['scalar_ell'],p['scalar_m'];a=args.spin;w=args.omega_cloud+(m-1)/(orbit**1.5+a)
   g=RadialGreen(a,p['alpha'],w,l,m,rmax=p['green_outer_radius'],offset=p['green_horizon_offset'],rtol=1e-11,infinity_method=p.get('infinity_method','series'))
   if g.infinity_method=='series' and g.series_last_term_relative>1e-3:raise ValueError('Unresolved boundary in fixed-source control')
   nodes=np.array([s['r'] for s in data['samples']]);J=np.array([complex(*s['source']) for s in data['samples']]);weights=np.array([s['weight'] for s in data['samples']]);oldrp=1+np.sqrt(1-p['metric']['a']**2)
   response=SampledResponse(OriginalSourceCoordinate(g,oldrp),p['source_panels'],nodes,J,log_first=p.get('horizon_log_first_panel',False))
   S=angular_mode(np.pi/2,l,m,a*a*(w*w-p['alpha']**2))[0];v=response.evaluate(r)[0][:,None]*S/p['alpha']**3*np.exp(1j*m*phi[None,:])
   values['all']+=v
   if (l,m) in POS:values['positive']+=v
   probes=np.array([nodes[0],3,20,orbit,100,320]);mode_details.append(dict(orbit=orbit,ell=l,m=m,omega=w,old_omega=p['omega'],wronskian_spread=float(np.max(abs(g.wronskian(probes)/g.w0-1)))))
  for label,field in values.items():
   for k,rr in enumerate(r):
    v=abs(field[k]);z=ref['nominal_abs_field'][io,k];lo=ref['lower_color_interval'][io,k];hi=ref['upper_color_interval'][io,k]
    rows.append(dict(orbit=orbit,selection=label,radius=float(rr),local_amplitude=v.tolist(),reference_amplitude=z.tolist(),rms_ratio=float(np.sqrt(np.mean(v*v)/np.mean(z*z))),normalized_amplitude_rms_difference=float(np.sqrt(np.mean((v-z)**2)/np.mean(z*z))),inside_color_interval_fraction=float(np.mean((v>=lo)&(v<=hi)))))
 result=dict(status='frozen_source_changed_propagator_control_not_new_background_solution',spin=args.spin,omega_cloud=args.omega_cloud,rows=rows,mode_details=mode_details,input_sha256=inputs,implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),limitations=['Metric source, background radial cloud and source angular projection remain at the original stationary spin and frequency.','Only all radial Green functions and reconstructed angular basis/frequency are replaced.','Original logarithmic interpolation coordinate is retained explicitly; no metric or cloud is evaluated outside its domain.','Reference is a raster color interval, not raw numerical data.'])
 args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 for row in rows:
  if row['selection']=='all' and row['radius'] in [50,100,150]:print({k:v for k,v in row.items() if k not in ['local_amplitude','reference_amplitude']})
if __name__=='__main__':main()
