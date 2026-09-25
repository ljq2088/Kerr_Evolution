"""Fast fixed-metric control isolating Leaver/BL-mass/order4 changes.

Uses ONLY the previously computed a=.88 scalar22 metrics at r0=41.1/42.1.
Low-L two-dimensional source is a diagnostic, not the final Li-style batch.
"""
import json,hashlib
from pathlib import Path
import numpy as np
from li_normalized_cloud import LiNormalizedCloud
from li_order4_green import LiOrder4Green
from environment_lorenz_mode import LorenzMetricMode
from environment_metric_sampling import PrecomputedMetric
from environment_source import project_source
from environment_cloud import mode_flux
from source_provenance import local_dependency_hashes
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    cloud=LiNormalizedCloud();rows=[]
    for rp,tag in [(41.1,'41p1'),(42.1,'42p1')]:
        p=ROOT/f'docs/root_cause_followup_20260917/fresh_a088_rp{tag}_sl2_sm2_L6_q12_nr8_h32.json'
        old=json.loads(p.read_text());par=old['parameters']
        base=LorenzMetricMode(rp,.88,1,6)
        assert par['metric']==base.provenance
        folder=ROOT/'outputs/metric_cache'/old['metric_sampling']['cache_key']
        theta=np.arccos(np.polynomial.legendre.leggauss(12)[0])
        radii=np.array([s['r'] for s in old['samples']]);weights=np.array([s['weight'] for s in old['samples']])
        values={};hashes={}
        for r in radii:
            f=folder/(hashlib.sha256(float(r).hex().encode()).hexdigest()[:24]+'.npz')
            with np.load(f,allow_pickle=False) as z:
                meta=json.loads(str(z['metadata']))
                assert meta['metric']==base.provenance and meta['source_hash']==old['metric_sampling']['source_hash']
                assert np.array_equal(np.array(meta['theta']),theta) and float(z['r'])==r
                for t,h in zip(theta,z['h']):values[(float(r),float(t))]=h
            hashes[f.name]=sha(f)
        metric=PrecomputedMetric(base,values)
        omega,J=project_source(cloud,radii,rp,2,2,metric,ntheta=12)
        oldJ=np.array([complex(*s['source']) for s in old['samples']])
        green=LiOrder4Green(.88,.3,omega,2,2,rmax=32000.)
        zi=np.sum(weights*green.insol.sol(radii)[0]*J)/green.w0
        zh=np.sum(weights*green.upsol.sol(radii)[0]*J)/green.w0
        flux=mode_flux(omega,2,cloud.omega,1,.3,.88,zi,zh)
        row=dict(orbit=rp,ell=2,m=2,baseline=str(p.relative_to(ROOT)),baseline_sha256=sha(p),
            metric_cache=folder.name,metric_cache_sha256=hashes,
            omega=omega,source_relative_l2=float(np.linalg.norm(J-oldJ)/np.linalg.norm(oldJ)),
            flux=flux,flux_Li_units={b:{k:v/.3**6 for k,v in F.items()} for b,F in flux.items()},
            old_flux=old['flux'],
            relative_change_percent={b:100*(flux[b]['orbital_energy']/old['flux'][b]['orbital_energy']-1)
                                     if old['flux'][b]['orbital_energy'] else None for b in flux},
            wronskian_relative_spread=float(np.max(abs(green.wronskian(radii)/green.w0-1))),
            infinity_boundary=green.boundary_audit)
        rows.append(row);print({k:v for k,v in row.items() if k not in ('metric_cache_sha256','infinity_boundary')},flush=True)
    out=ROOT/'docs/li_alignment/fixed_metric_flux_control.json'
    out.write_text(json.dumps(dict(status='completed_diagnostic_not_full_Li_comparison',rows=rows,
        cloud=cloud.provenance,mass_audit=cloud.normalization_audit(),
        implementation_sha256=local_dependency_hashes(ROOT/'src',['report_li_fixed_metric_flux_control']),
        limitations=['Only scalar22; no total flux or direct comparison to a total curve.',
                     'Existing fresh-a=.88 L6 q12 metric; source is direct two-dimensional contraction.',
                     'Main j18 separated-source r0=20 batch is independent of this diagnostic.']),indent=2)+'\n')
if __name__=='__main__':main()
