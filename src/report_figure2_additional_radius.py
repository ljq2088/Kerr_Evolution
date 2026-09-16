"""Publish a provenance-checked partial Figure-2 radius without calling it total."""
import argparse,hashlib,json
import numpy as np
from pathlib import Path
from source_provenance import source_fingerprint,validate_saved_samples
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'docs/environment_reproduction'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--radius',type=float,default=15.);p.add_argument('--reference',type=Path,required=True);p.add_argument('--published-reference',type=Path);p.add_argument('--coverage',type=Path,required=True);p.add_argument('--batch',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
 ref=json.loads(args.reference.read_text());coverage=json.loads(args.coverage.read_text());batch=json.loads(args.batch.read_text());r=args.radius
 if ref['status']!='digitized_reference_not_author_numerical_data' or ref['background']!='kerr':raise ValueError('Require original Kerr vector-PDF extraction')
 target=next(row for row in ref['rows'] if row['rp']==r)
 if batch['status']!='batch_completed_finite_resolution_not_converged':raise ValueError('Batch not finished')
 if coverage['parameters']['rp']!=r or coverage['parameters']['metric_ellmax']!=18:raise ValueError('Wrong coverage parameters')
 actuals={};input_hashes={str(p.relative_to(ROOT)):sha(p) for p in [args.reference.resolve(),args.coverage.resolve(),args.batch.resolve()]};channels=[];fingerprint=source_fingerprint()
 for entry in batch['channels']:
  path=FOLDER/entry['file'];data=json.loads(path.read_text());validate_saved_samples(data,fingerprint)
  par=data['parameters'];key=(par['scalar_ell'],par['scalar_m'])
  if par['metric']['orbital_radius']!=r or par['metric']['ellmax']!=18 or par['alpha']!=.3 or par.get('background') is not None:raise ValueError('Channel physics mismatch')
  if key in actuals:raise ValueError('Repeated channel')
  actuals[key]=data;input_hashes[str(path.relative_to(ROOT))]=sha(path)
  channels.append(dict(ell=key[0],m=key[1],metric_m=par['metric']['m_g'],metric_conjugated=entry['metric_conjugated'],file=entry['file'],status=data['status'],sample_count=len(data['samples']),source_provenance_sha256=data['source_provenance']['sha256'],z_h=data['z_h'],z_inf=data['z_inf'],flux=data['flux'],wronskian_relative_spread=data['wronskian_relative_spread'],infinity_boundary_audit=data.get('infinity_boundary_audit'),outer_source_cutoff_sequence=data.get('outer_source_cutoff_sequence')))
 published=None
 if args.published_reference is not None:
  published=json.loads(args.published_reference.read_text())
  if published['status']!='digitized_reference_only_not_local_simulation' or published['normalization']!='per q^2 epsilon^2':raise ValueError('Unexpected published reference convention')
  input_hashes[str(args.published_reference.resolve().relative_to(ROOT))]=sha(args.published_reference)
 boundaries={}
 for name in ('infinity','horizon'):
  present=[];missing=[];values=[]
  for item in coverage[name]['modes']:
   key=(item['ell'],item['m'])
   if 'flux' in item:
    if key not in actuals or item['flux']!=actuals[key]['flux']:raise ValueError('Coverage includes unexpected or stale flux')
    present.append(list(key));values.append(actuals[key]['flux'][name]['orbital_energy'])
   else:missing.append(list(key))
  if len(present)!=coverage[name]['computed_count']:raise ValueError('Coverage count mismatch')
  partial=sum(values);reference=target['flux']['infinity' if name=='infinity' else 'horizon_magnitude']['per_q2_cloud_mass'];complete=not missing
  if not complete and coverage[name]['finite_resolution_total'] is not None:raise ValueError('Incomplete coverage improperly has total')
  boundaries[name]=dict(required_count=coverage[name]['required_count'],computed_count=len(present),computed_modes=present,missing_modes=missing,partial_signed_flux=partial,total_signed_flux=partial if complete else None,is_total=complete,historical_v1_total_magnitude=reference,partial_magnitude_divided_by_historical_v1_total=abs(partial)/reference,comparison_kind='partial sum versus explicitly versioned reference total; ratio is not a total-flux accuracy estimate',converged=False)
  comparisons={'dyson_arxiv_2501_09806v1':dict(total_magnitude=reference,partial_magnitude_divided_by_reference_total=abs(partial)/reference,reference_file=str(args.reference),reference_role='historical arXiv v1; not the later corrected curve')}
  if published is not None:
   for prefix,label in [('dyson','later_dyson_curve_in_li_v2'),('li','li_2507_02045v2')]:
    curve=published['curves'][prefix+'_'+name];xx=np.asarray([v['rp_over_M'] for v in curve]);yy=np.asarray([v['flux_magnitude'] for v in curve])
    if not np.all(np.diff(xx)>0) or not xx[0]<=r<=xx[-1] or not np.all(yy>0):raise ValueError('Invalid reference curve or unsupported radius')
    mag=float(np.exp(np.interp(r,xx,np.log(yy)))*.3**6)
    comparisons[label]=dict(total_magnitude=mag,partial_magnitude_divided_by_reference_total=abs(partial)/mag,reference_file=str(args.published_reference),interpolation='linear in log flux vs radius; converted by alpha^6 to per q^2 (Mc/M)')
  boundaries[name]['versioned_reference_comparisons']=comparisons
 result=dict(status='fresh_additional_radius_partial_modes_not_paper_total_reproduction',rp=r,alpha=.3,normalization='per q^2 (Mc/M); no fitted rescaling',parameters=coverage['parameters'],boundaries=boundaries,channels=channels,metric_precomputation=batch.get('metric_precomputation'),input_sha256=input_hashes,source_provenance=fingerprint,implementation_sha256=sha(Path(__file__)),limitations=['Only the explicitly listed metric m pair and scalar multipoles were calculated. Missing modes are not set to zero.','New radius uses fresh current-source metric data; no historical tensor cache is relabelled after source-hash changes.','Finite Lg18/nr8/nt18 data are not a complete convergence demonstration.','Versioned references are vector-PDF interpolations, not author numerical tables, with no formal digitization error bound.','The old arXiv v1 horizon curve is retained as historical evidence and must not be conflated with the later Dyson curve or Li curve.','This is not a new total Figure-2 point; plotting requires a distinct partial marker.'])
 cache_audits={}
 if r==15.:
  for name in ('figure2_rp15_auxiliary_cache_20260916.json','figure2_rp15_cache_handoff_20260916.json'):
   path=FOLDER/name;audit=json.loads(path.read_text());source_hash=audit.get('source_hash',audit.get('metadata',{}).get('source_hash'))
   if source_hash!=result['metric_precomputation']['source_hash']:raise ValueError('Cache generation source hash differs from final batch')
   cache_audits[name]=dict(status=audit['status'],sha256=sha(path),source_hash=source_hash);result['input_sha256'][str(path.relative_to(ROOT))]=sha(path)
 result['cache_generation_audits']=cache_audits
 args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(boundaries,indent=2))
if __name__=='__main__':main()
