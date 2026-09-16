"""Assemble provenance and output locations for the two fresh Li static sectors."""
import hashlib,json
from pathlib import Path
from source_provenance import source_fingerprint,validate_saved_samples
ROOT=Path(__file__).resolve().parents[1];FOLDER=ROOT/'docs/environment_reproduction'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 current=source_fingerprint();rows=[]
 for rp in (41.1,42.1):
  sample=ROOT/f'outputs/static_samples_a0.877153_r{rp:g}_q12_eps5e-05.npz'
  match=FOLDER/f'static_tetrad_a0.877153_r{rp:g}_L6_q12_j8_free8_paper.json';md=json.loads(match.read_text());match_argument=str(match.relative_to(ROOT))
  candidates=[]
  for path in FOLDER.glob('scalar_batch_mg0_*.json'):
   data=json.loads(path.read_text());p=data.get('parameters',{})
   if p.get('orbital_radius')==rp and p.get('static_matching')==match_argument and p.get('metric_ellmax')==6 and p.get('radial_order')==8 and p.get('angular_order')==12 and p.get('scalar_ells')==[3,5]:candidates.append((path,data))
  if len(candidates)!=1:raise ValueError(f'Expected one unique static batch at {rp}: {len(candidates)}')
  path,batch=candidates[0];channels=[]
  for item in batch.get('channels',[]):
   response=FOLDER/item['file'];rd=json.loads(response.read_text());validate_saved_samples(rd,current);par=rd['parameters'];samples=rd['samples']
   if par['scalar_m']!=1 or par['scalar_ell'] not in (3,5) or par['metric']['orbital_radius']!=rp or par['metric']['ellmax']!=6:raise ValueError('Wrong static scalar data')
   orbit=next(v for v in rd['radial_response_at_panel_boundaries'] if v['r']==rp)
   channels.append(dict(file=response.name,sha256=sha(response),ell=par['scalar_ell'],m=par['scalar_m'],status=rd['status'],source_samples=len(samples),source_provenance_sha256=rd['source_provenance']['sha256'],radial_at_particle=orbit['field'],radial_derivative_at_particle=orbit['derivative'],wronskian_relative_spread=rd['wronskian_relative_spread'],flux=rd['flux']))
  row=dict(rp=rp,status=batch['status'],batch_file=path.name,batch_sha256=sha(path),sample_file=str(sample.relative_to(ROOT)),sample_sha256=sha(sample),matching_file=match.name,matching_sha256=sha(match),matching_a=md['sample_metadata']['a'],matching_rank=md['rank'],matching_selected_continuity_residual=md['maximum_selected_continuity_residual'],matching_value_residual_by_degree=md['scaled_value_residual_by_degree'],matching_derivative_residual_by_degree=md['scaled_derivative_residual_by_degree'],metric_precomputation=batch.get('metric_precomputation'),channels=channels)
  rows.append(row)
 complete=all(r['status']=='batch_completed_finite_resolution_not_converged' and len(r['channels'])==2 for r in rows)
 result=dict(status='four_static_channels_completed_for_li_orbits_finite_L6' if complete else 'static_channels_in_progress',current_source_provenance=current,parameters=dict(alpha=.3,a=.8771530275949366,metric_m=0,scalar_m=1,scalar_ells=[3,5],metric_ellmax=6,angular_order=12,radial_order=8,horizon_order=32,source_inner_offset=.0005,source_outer_radius=320.,green_outer_radius=4000.,infinity_method='coulomb',workers_per_orbit=3),rows=rows,input_implementation_sha256={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'src/report_static_matching.py',ROOT/'src/report_static_tetrad_matching.py',ROOT/'src/report_environment_mode_batch.py',Path(__file__)]},limitations=['Static matching uses paper independent continuity constraints and exact conserved charges; high-degree truncation-edge derivative residuals remain explicitly recorded.','These two static m=1 channels per orbit complete only the static subset of the intended ell2..5 field comparison.','Finite L6/qtheta12/nr8 output is not a complete convergence demonstration.','Threshold static channels affect the complex field despite zero orbital-effective scalar flux.'])
 out=FOLDER/'li_new_orbits_static_20260917.json';out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(status=result['status'],rows=[dict(rp=r['rp'],status=r['status'],channels=[c['file'] for c in r['channels']]) for r in rows]),indent=2))
if __name__=='__main__':main()
