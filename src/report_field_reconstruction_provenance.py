"""Audit provenance and available like-for-like refreshed spatial-field inputs."""
import hashlib,json
from pathlib import Path
from collections import Counter
import numpy as np
from source_provenance import source_fingerprint,samples_hash
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'docs/environment_reproduction'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def identity(data):
    p=data['parameters'];m=p['metric'];return (p['alpha'],m['a'],m['orbital_radius'],p['cloud_mass'],p['scalar_ell'],p['scalar_m'],m['m_g'])
def z(a):return complex(*a)
def main():
    current=source_fingerprint();candidate_index={}
    for path in D.glob('*.json'):
        raw=path.read_bytes()
        if b'"source_provenance"' not in raw:continue
        d=json.loads(raw)
        if d.get('status')!='truncated_single_mode_not_converged' or not isinstance(d.get('source_provenance'),dict) or d.get('flux_validity'):continue
        if d.get('source_samples_sha256')!=samples_hash(d['samples']):continue
        candidate_index.setdefault(identity(d),[]).append((path,d))
    reports=[]
    for name in ['figure5_threshold_wakes.json','wake_diagnostic_rp20_L18_sl12_mixed_wide_linear_eps.json','wake_diagnostic_rp3.5_L18_sl12_mixed_wide_linear_eps.json']:
        path=D/name;fig=json.loads(path.read_text())
        for ds in fig.get('datasets',[fig]):
            counts=Counter();inputs=[]
            for item in ds['inputs']:
                file=Path(item['file']);file=D/file if len(file.parts)==1 else ROOT/file
                d=json.loads(file.read_text());counts['inputs']+=1;versioned=isinstance(d.get('source_provenance'),dict);counts['versioned' if versioned else 'unversioned']+=1;hash_matches=sha(file)==item['sha256'];counts['saved_input_hash_matches']+=hash_matches
                row=dict(file=str(file.relative_to(ROOT)),sha256=sha(file),saved_plot_input_hash_matches=hash_matches,ell=d['parameters']['scalar_ell'],m=d['parameters']['scalar_m'],versioned=versioned,refresh_candidates=[])
                for newfile,new in candidate_index.get(identity(d),[]):
                    if newfile==file:continue
                    changes={k:[v,new['parameters'].get(k)] for k,v in d['parameters'].items() if new['parameters'].get(k)!=v};oldpt=next(s for s in d['radial_response_at_panel_boundaries'] if s['r']==ds['parameters']['r0']);newpt=next(s for s in new['radial_response_at_panel_boundaries'] if s['r']==ds['parameters']['r0'])
                    changedfiles=[k for k,v in new['source_provenance']['files'].items() if current['files'].get(k)!=v]
                    compare=dict(file=str(newfile.relative_to(ROOT)),sha256=sha(newfile),same_discretization=(not changes),same_source_discretization=not(set(changes)-{'green_outer_radius','green_horizon_offset','infinity_method'}),parameter_changes=changes,source_samples_hash_valid=True,current_source_file_mismatches=changedfiles,current_runtime_matches=new['source_provenance']['runtime']==current['runtime'],version_scope='current' if new['source_provenance']==current else 'versioned historical; original hashes retained',relative_orbit_complex_field_change=abs(z(newpt['field'])-z(oldpt['field']))/max(abs(z(oldpt['field'])),1e-300),relative_ZH_change=abs(z(new['z_h'])-z(d['z_h']))/max(abs(z(d['z_h'])),1e-300))
                    if compare['same_source_discretization'] and len(d['samples'])==len(new['samples']):
                        assert [s['r'] for s in d['samples']]==[s['r'] for s in new['samples']]
                        j0=np.array([z(s['source']) for s in d['samples']]);j1=np.array([z(s['source']) for s in new['samples']]);compare['relative_source_max_norm']=float(np.max(abs(j1-j0))/np.max(abs(j0)))
                    row['refresh_candidates'].append(compare)
                if row['refresh_candidates']:counts['inputs_with_versioned_refresh']+=1
                inputs.append(row)
            reports.append(dict(plot=str(path.relative_to(ROOT)),plot_sha256=sha(path),orbit=ds['parameters']['r0'],counts=dict(counts),inputs=inputs))
    # The notorious threshold comparison changes several controls simultaneously.
    names=['forced_mode_nr8_nt10_L6_alpha0.3_rp41.6_mg1_sl2_gh0.0001_go32000_coulomb_log_h16.json','forced_mode_nr8_nt18_L18_alpha0.3_rp41.6_mg1_sl2_gh0.0001_go32000_coulomb_inner0.0005_outer320_log_h32.json'];a,b=[json.loads((D/n).read_text()) for n in names]
    changes={k:[v,b['parameters'].get(k)] for k,v in a['parameters'].items() if b['parameters'].get(k)!=v};paired=[]
    for name in ['forced_mode_nr8_nt18_L18_alpha0.3_rp41.6_mg1_sl2_gh0.0001_go32000_coulomb_inner0.0005_outer640_log_h32.json','forced_mode_nr8_nt18_L18_alpha0.3_rp41.6_mg1_sl2_gh0.0001_go64000_coulomb_inner0.0005_outer320_log_h32.json']:
        v=json.loads((D/name).read_text());paired.append(dict(file=name,sha256=sha(D/name),relative_FH_change=v['flux']['horizon']['orbital_energy']/b['flux']['horizon']['orbital_energy']-1,relative_ZH_change=abs(z(v['z_h'])-z(b['z_h']))/abs(z(b['z_h'])),parameter_changes={k:[x,v['parameters'].get(k)] for k,x in b['parameters'].items() if v['parameters'].get(k)!=x},provenance='unversioned historical'))
    result=dict(status='field_input_version_audit_not_physical_source_validation',plots=reports,threshold_L6_to_L18=dict(inputs=[dict(file=n,sha256=sha(D/n)) for n in names],parameter_changes=changes,FH=[a['flux']['horizon']['orbital_energy'],b['flux']['horizon']['orbital_energy']],relative_FH_change=b['flux']['horizon']['orbital_energy']/a['flux']['horizon']['orbital_energy']-1,interpretation='Source inner offset0.05->0.0005, horizon Gauss16->32, angular10->18 and metricL6->L18 changed together; not a pure metric-cutoff or source-version effect.'),threshold_same_L18_controls=paired,conclusions=['All four examined plot datasets consist of88 unversioned historical source responses each; file hashes bind saved plots to exact inputs but cannot prove the source-code version.','At rp20, five spatial-field inputs have a versioned fresh recomputation with the same source grid from the previous audit; Green outer radius is1000 instead of the plotting input4000, and this difference is recorded explicitly. Their pre-guard source fingerprints are retained, not relabeled as current.','No versioned rp41.6/rp41.8 orrp3.5 field refresh is present; independent Green reconstruction alone cannot certify those original source samples.','Historical threshold radius/boundary refinements may inform sensitivity but are not promoted to present-code provenance.'],implementation_sha256=sha(Path(__file__)))
    p=D/'field_reconstruction_provenance_20260916.json';p.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for r in reports:
        print(r['orbit'],r['counts'])
        for v in r['inputs']:
            for c in v['refresh_candidates']:print('refresh',v['ell'],v['m'],c['relative_orbit_complex_field_change'],c['same_discretization'])
    print('threshold joint changes',result['threshold_L6_to_L18']['relative_FH_change'])
if __name__=='__main__':main()
