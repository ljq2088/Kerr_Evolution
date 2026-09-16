"""Inventory and quantify trustworthy partial refreshes of historical Fig.7 data.

No response is regenerated, rehashed as a new calculation, or silently promoted.
"""
import hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from environment_source import angular_mode
from source_provenance import source_fingerprint,validate_saved_samples
ROOT=Path(__file__).resolve().parents[1];FOLDER=ROOT/'docs/environment_reproduction'
GUARD='    # OdeSolution otherwise silently extrapolates beyond the integrated domain.\n    if not np.isfinite(r) or not radial.rmin <= r <= radial.rmax:\n        raise ValueError(f"Trace/kappa requested r={r!r} outside solved radial range "\n                         f"[{radial.rmin!r}, {radial.rmax!r}]; radius must be finite")\n'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pair(z):return [float(z.real),float(z.imag)]
def main():
 current=source_fingerprint();kappa=(ROOT/'src/lorenz_kappa.py').read_text()
 # Guard text is fixed explicitly, not selected through a permissive regex.
 guard=GUARD.replace('\\n','\n')
 if kappa.count(guard)!=1:raise ValueError('Expected precise domain guard once')
 old_kappa_hash=hashlib.sha256(kappa.replace(guard,'').encode()).hexdigest()
 grouped=defaultdict(list);records=[]
 for path in sorted(FOLDER.glob('*.json')):
  try:data=json.loads(path.read_text())
  except (ValueError,OSError):continue
  if not isinstance(data,dict):continue
  p=data.get('parameters',{})
  if not isinstance(p,dict):continue
  g=p.get('metric',{})
  if not isinstance(g,dict) or p.get('alpha')!=.3 or g.get('orbital_radius')!=20 or g.get('a')!=.8771530275949366 or 'scalar_ell' not in p or 'radial_response_at_panel_boundaries' not in data:continue
  points=[r for r in data['radial_response_at_panel_boundaries'] if r['r']==20]
  if len(points)!=1:continue
  prov=data.get('source_provenance');classification='unversioned_historical';checksum_ok=False;bounds=None
  if prov:
   try:validate_saved_samples(data,prov);checksum_ok=True
   except ValueError:classification='source_sample_checksum_failure'
   if checksum_ok:
    if prov==current:classification='exact_current_source'
    else:
     rewritten=json.loads(json.dumps(prov));rewritten.pop('sha256',None);expected=json.loads(json.dumps(current));expected.pop('sha256',None)
     if rewritten['files'].get('lorenz_kappa')==old_kappa_hash:
      rewritten['files']['lorenz_kappa']=current['files']['lorenz_kappa']
      radii=np.array([z['r'] for z in data['samples']]);rplus=1+np.sqrt(1-g['a']**2)
      bounds=dict(min_source_r=float(radii.min()),max_source_r=float(radii.max()),kappa_rmin=float(rplus+1e-4),kappa_rmax=2000.)
      if rewritten==expected and np.all(np.isfinite(radii)) and np.all((radii>=rplus+1e-4)&(radii<=2000.)):classification='pre_domain_guard_same_numerical_algorithm_in_domain'
      else:classification='different_versioned_source'
     else:classification='different_versioned_source'
  row=dict(file=path.name,sha256=sha(path),ell=p['scalar_ell'],m=p['scalar_m'],metric_ellmax=g.get('ellmax'),radial_order=p.get('radial_order'),angular_order=p.get('angular_order'),green_outer_radius=p.get('green_outer_radius'),status=data.get('status'),flux_validity=data.get('flux_validity'),sample_count=len(data.get('samples',[])),source_sha256=prov.get('sha256') if prov else None,classification=classification,source_sample_checksum_valid=checksum_ok,kappa_guard_domain_check=bounds,radial_at_particle=points[0]['field'],mtime=path.stat().st_mtime)
  records.append(row);grouped[(row['ell'],row['m'])].append((row,data))
 particle_path=FOLDER/'particle_field_L18.json';particle=json.loads(particle_path.read_text());markers={v['ell']:v['plotted_value'] for v in json.loads((FOLDER/'paper_figure7_markers.json').read_text())['markers']};shells=[]
 for shell in particle['multipoles']:
  if not shell['complete']:continue
  ell=shell['ell'];before=complex(*shell['particle_field_sum']);after=before;used=[];missing=[];exact=[]
  for component in shell['components']:
   m=component['m'];eligible=[(r,d) for r,d in grouped[(ell,m)] if r['classification'] in ('exact_current_source','pre_domain_guard_same_numerical_algorithm_in_domain') and r['metric_ellmax']==18 and r['status']=='truncated_single_mode_not_converged' and r['flux_validity']!='historical_unreliable_boundary_result']
   if not eligible:missing.append(m);continue
   row,new=max(eligible,key=lambda v:(v[0]['classification']=='exact_current_source',v[0]['angular_order'],v[0]['radial_order'],v[0]['mtime']))
   oldpath=FOLDER/component['file'];old=json.loads(oldpath.read_text());p=new['parameters'];alpha=p['alpha']
   s=float(angular_mode(np.pi/2,ell,m,p['metric']['a']**2*(p['omega']**2-alpha**2))[0]);zold=complex(*component['radial']);z=complex(*row['radial_at_particle']);delta=(z-zold)*s;after+=delta
   node_test=dict(same_radial_nodes=False)
   if len(old['samples'])==len(new['samples']) and all(a['r']==b['r'] for a,b in zip(old['samples'],new['samples'])):
    a=np.array([complex(*v['source']) for v in old['samples']]);b=np.array([complex(*v['source']) for v in new['samples']]);node_test=dict(same_radial_nodes=True,source_relative_L2=float(np.linalg.norm(b-a)/np.linalg.norm(a)),source_max_absolute=float(np.max(abs(b-a))))
   used.append(dict(m=m,old_file=component['file'],old_sha256=sha(oldpath),new_file=row['file'],new_sha256=row['sha256'],classification=row['classification'],radial_relative_change=float(abs(z-zold)/abs(zold)),old_green_outer_radius=old['parameters'].get('green_outer_radius'),new_green_outer_radius=p.get('green_outer_radius'),angular_weight=s,particle_field_delta=pair(delta),source_comparison=node_test))
   if row['classification']=='exact_current_source':exact.append(m)
  shells.append(dict(ell=ell,required_m=[v['m'] for v in shell['components']],trusted_m=[v['m'] for v in used],exact_current_m=exact,missing_trusted_m=missing,complete_trusted_shell=not missing,historical_field=pair(before),mixed_partial_refresh_field=pair(after),complex_relative_change=float(abs(after-before)/abs(before)),magnitude_fractional_change=float(abs(after)/abs(before)-1),historical_over_paper_marker=float(abs(before)/.3**3/markers[ell]),mixed_partial_refresh_over_paper_marker=float(abs(after)/.3**3/markers[ell]),replacement_modes=used))
 result=dict(status='no_complete_current_or_guard_equivalent_fig7_shell_available_partial_refresh_quantified',current_source_provenance=current,guard_equivalence_proof=dict(current_kappa_sha256=sha(ROOT/'src/lorenz_kappa.py'),guard_removed_sha256=old_kappa_hash,exact_guard_text=guard,method='Delete exactly this guard from current file and compare hash; all other fingerprint fields must match; require all source radii inside unchanged default trace/kappa radial interval. Original input provenance is never changed.'),inventory_counts=dict(Counter(r['classification'] for r in records)),records=records,particle_input_sha256=sha(particle_path),shells=shells,summary=dict(required_particle_channels=sum(len(s['required_m']) for s in shells),trusted_guard_equivalent_or_current_channels=sum(len(s['trusted_m']) for s in shells),exact_current_channels=sum(len(s['exact_current_m']) for s in shells),complete_trusted_shells=[s['ell'] for s in shells if s['complete_trusted_shell']],max_partial_refresh_complex_relative_change=max(s['complex_relative_change'] for s in shells)),implementation_sha256=sha(Path(__file__)),limitations=['Partial replacements retain all other historical channels; these mixed sums are a sensitivity diagnostic, not updated certified Fig.7 points.','Pre-guard data keep their original source hash; byte-level guard-only equivalence is not source-hash equality.','Different Green outer boundaries can also contribute to the small changes; this is not a clean one-parameter algorithm comparison.','Missing trusted modes require actual new calculations; lack of provenance alone is not evidence that their numbers are wrong.'])
 out=FOLDER/'particle_field_source_inventory_20260916.json';out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result['summary'],indent=2));print(json.dumps(result['inventory_counts'],indent=2));print(json.dumps([{k:s[k] for k in ['ell','trusted_m','missing_trusted_m','complex_relative_change','magnitude_fractional_change']} for s in shells],indent=2))
if __name__=='__main__':main()
