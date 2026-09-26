"""Fresh Dyson Fig.1 geometry with Li-aligned fixed-spin cloud/source conventions.

No polling service: this process computes, checkpoints and renders sequentially.
Figures are finite-truncation calculations, not claims of matching author arrays.
"""
import argparse,json,os,subprocess,sys,traceback
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from li_normalized_cloud import LiNormalizedCloud
from li_field_green import LiFieldGreen
from li_separated_source import LiSeparatedSource
from environment_dense_metric import DenseLorenzMetricMode
from environment_static_lorenz import StaticLorenzMode
from environment_metric_sampling import precompute_metric
from environment_response import SampledResponse
from environment_source import angular_mode
from environment_cloud import mode_flux
from report_environment_forced_mode import source_grid
from report_li_aligned_flux import spherical_coefficients,enc,sha,stamp,save
from source_provenance import local_dependency_hashes
ROOT=Path(__file__).resolve().parents[1]

def field_inventory(lmax=12):
    return [(l,m) for l in range(2,lmax+1) for m in range(-l,l+1) if (l+m)%2==0]

def field_grid(cloud,r0,mg,nr=24,nh=64):
    cfg=SimpleNamespace(orbital_radius=r0,source_inner_offset=5e-4,
        source_outer_radius=320.,radial_order=nr,horizon_order=nh,horizon_log=True)
    panels,_,_=source_grid(cfg,cloud)
    wg=mg/(r0**1.5+cloud.a);w=cloud.omega+wg
    k=np.sqrt(max(w*w-cloud.mu**2,0.))
    # Resolve oscillatory radial products; at most six phase cycles per panel.
    width=min(40.,12*np.pi/max(abs(wg)+k,.01))
    refined=[panels[0],panels[1]]
    for lo,hi in zip(panels[1:-1],panels[2:]):
        refined.extend(np.linspace(lo,hi,int(np.ceil((hi-lo)/width))+1)[1:])
    panels=np.asarray(refined);nodes=[];weights=[]
    for i,(lo,hi) in enumerate(zip(panels[:-1],panels[1:])):
        x,w=np.polynomial.legendre.leggauss(nh if i==0 else nr)
        if i==0:
            l,h=np.log(lo-cloud.rp),np.log(hi-cloud.rp)
            d=np.exp((l+h)/2+(h-l)*x/2);r=cloud.rp+d;wt=(h-l)/2*w*d
        else:r=(lo+hi)/2+(hi-lo)*x/2;wt=(hi-lo)/2*w
        nodes.extend(r);weights.extend(wt)
    return panels,np.array(nodes),np.array(weights)

def static_matching(out):
    dest=out/'static_matching.json'
    if dest.exists():return json.loads(dest.read_text())
    command=[sys.executable,'-u','src/report_static_matching.py','--a','.88','--r0','3.5',
        '--ellmax','20','--quadrature','48','--testmax','20','--continuity-only','--reuse-samples']
    with (out/'static_matching.log').open('a') as stream:
        subprocess.run(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,check=True)
    p=ROOT/'docs/environment_reproduction/static_matching_a0.88_r3.5_L20_q48_j20_continuity_exactcharges.json'
    z=ROOT/'outputs/static_samples_a0.88_r3.5_q48_eps5e-05.npz'
    d=json.loads(p.read_text())
    with np.load(z,allow_pickle=False) as bank:d['sample_metadata']=json.loads(str(bank['metadata']))
    d['raw_matching_sha256']=sha(p);d['sample_bank_sha256']=sha(z)
    save(dest,d);return d

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);ap.add_argument('--workers',type=int,default=4)
    args=ap.parse_args()
    if not 1<=args.workers<=4:raise ValueError('Use 1..4 workers')
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    manifest=out/'execution.json'
    fingerprints=local_dependency_hashes(ROOT/'src',['report_figure1_background','render_figure1_background'])
    fingerprints['li_source_terms.json']=sha(ROOT/'src/li_source_terms.json')
    fingerprint0=dict(fingerprints)
    cloud=LiNormalizedCloud();expected=field_inventory()
    config=dict(target='Dyson arXiv:2501.09806v1 Fig.1',orbit=3.5,a=.88,mu=.3,
        cloud=cloud.provenance,scalar_lmin=2,scalar_lmax=12,expected_modes=expected,
        metric_input_L=20,metric_output_j=18,metric_q=40,source_q=64,pmax=12,
        radial_order=24,horizon_order=64,source_outer=320.,source_inner_offset=5e-4,
        static_completion='berndtson',field_scale='unit cloud mass amplitude / alpha^3',
        positive_negative_static_and_bound_modes=True,
        coordinate_map='BL-label x=r sin(theta)cos(phi),y=r sin(theta)sin(phi),z=r cos(theta); paper map unverified',
        time=0,amplitude_or_phase_fitted=False)
    # Normalize tuples for exact JSON resume comparisons.
    config=json.loads(json.dumps(config))
    if manifest.exists():
        old=json.loads(manifest.read_text())
        if old['config']!=config or old['implementation_sha256']!=fingerprints:
            raise ValueError('Refuse incompatible resume')
    record=dict(status='running',pid=os.getpid(),started_utc=stamp(),config=config,
        implementation_sha256=fingerprints,completed_modes=[],expected_mode_count=len(expected),
        limitations=['Finite resolution; all retained modes do not imply convergence.',
         'Dyson caption a=.88 versus synchronized background in its method: use exact .88 consistently here.',
         'Author complex-frequency policy, coordinate embedding and static completion require verification.',
         'Local metric implementation, not original author data at this parameter point.'])
    def persist():
        record['updated_utc']=stamp();save(manifest,record)
    persist()
    try:
        x,wt=np.polynomial.legendre.leggauss(40);theta=np.arccos(x)
        rfield=np.unique(np.r_[np.geomspace(cloud.rp+.05,20.,181),np.linspace(20.,192.,1501)])
        for absmg in range(14):
            targets=[(l,m) for l,m in expected if abs(m-1)==absmg]
            pending=[]
            for l,m in targets:
                p=out/f'mode_l{l}_m{m}.json'
                if p.exists():
                    d=json.loads(p.read_text())
                    if d['implementation_sha256']!=fingerprint0 or d['config']!=config:
                        raise ValueError('Saved mode differs')
                    if sha(out/d['radial_file'])!=d['radial_sha256']:raise ValueError('Radial checksum mismatch')
                    record['completed_modes'].append([l,m])
                else:pending.append((l,m))
            persist()
            if not pending:continue
            record['stage']='static_matching' if absmg==0 else 'metric_sampling'
            record['metric_abs_m']=absmg;persist()
            base=StaticLorenzMode(static_matching(out)) if absmg==0 else DenseLorenzMetricMode(3.5,.88,absmg,20)
            # Common grid for the two metric-conjugate sources, using tighter spacing.
            plus=field_grid(cloud,3.5,absmg);minus=field_grid(cloud,3.5,-absmg)
            panels,radii,weights=plus if len(plus[1])>=len(minus[1]) else minus
            print(f'METRIC mg={absmg}, radii={len(radii)}',flush=True)
            metric,audit=precompute_metric(base,radii,theta,ROOT/'outputs/metric_cache',args.workers)
            now=local_dependency_hashes(ROOT/'src',['report_figure1_background','render_figure1_background'])
            now['li_source_terms.json']=sha(ROOT/'src/li_source_terms.json')
            if now!=fingerprint0:raise ValueError('Implementation changed while metric workers ran')
            record['stage']='source_and_field';persist()
            for mg in sorted({m-1 for _,m in pending}):
                bank=np.array([spherical_coefficients(metric,r,theta,wt,.88,mg,18) for r in radii])
                bp=out/f'metric_mg{mg}.npz'
                np.savez_compressed(bp,radii=radii,coefficients=bank,metadata=json.dumps(audit))
                for ell,m in [v for v in pending if v[1]-1==mg]:
                    omega=cloud.omega+mg/(3.5**1.5+.88)
                    target=lambda t:angular_mode(t,ell,m,.88**2*(omega**2-.3**2))[0]
                    projector=LiSeparatedSource(a=.88,omega_c=cloud.omega,m_c=1,metric_m=mg,
                        spherical_lmax=18,cloud_angular=cloud.angular_state,target_angular=target,
                        quadrature=64,pmax=12,dps=64)
                    J=np.array([projector.project(r,cloud.radial_state(r)[:,0],c) for r,c in zip(radii,bank)])
                    gp=LiFieldGreen(.88,.3,omega,ell,m)
                    spread=float(np.max(abs(gp.wronskian(radii)/gp.w0-1)))
                    if not np.isfinite(spread) or spread>1e-5:raise ValueError(f'Wronskian failure {ell,m}: {spread}')
                    response=SampledResponse(gp,panels,radii,J,log_first=True)
                    values=response.evaluate(rfield)
                    rp=out/f'radial_l{ell}_m{m}.npz'
                    np.savez_compressed(rp,r=rfield,field=values[0],derivative=values[1])
                    row=dict(status='finite_mode_complete',ell=ell,m=m,omega=omega,config=config,
                        implementation_sha256=fingerprint0,source_radii=radii.tolist(),
                        source_weights=weights.tolist(),source_panels=panels.tolist(),source=enc(J),
                        metric_bank_sha256=sha(bp),metric_sampling=audit,
                        radial_file=rp.name,radial_sha256=sha(rp),
                        wronskian_relative_spread=spread,green_outer=gp.rmax,boundary=gp.boundary_audit,
                        flux=mode_flux(omega,m,cloud.omega,1,.3,.88,response.up_coefficient,response.horizon_coefficient),
                        completed_utc=stamp())
                    save(out/f'mode_l{ell}_m{m}.json',row)
                    record['completed_modes'].append([ell,m]);persist()
                    print(f'FIELD {ell,m}: {len(record["completed_modes"])}/88',flush=True)
        if set(map(tuple,record['completed_modes']))!=set(expected):raise ValueError('Missing field modes')
        record['stage']='rendering_and_reference_comparison';persist()
        subprocess.run([sys.executable,'src/render_figure1_background.py',str(out)],cwd=ROOT,check=True)
        record.update(status='completed_finite_resolution_field_convergence_pending',completed_utc=stamp())
        record.pop('stage',None);persist()
    except Exception as exc:
        record.update(status='failed_resumable',exception=repr(exc),traceback=traceback.format_exc());persist();raise
if __name__=='__main__':main()
