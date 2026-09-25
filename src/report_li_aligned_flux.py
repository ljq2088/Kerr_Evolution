"""Resumable fresh a=.88 Li-style finite-resolution flux calculation.

Full output spherical j<=18 is projected from independently reconstructed
Lorenz spheroidal input L<=20. This is not an original-author code run.
All totals carry completeness checks and explicit author/convergence gaps.
"""
import argparse,hashlib,json,os,time,traceback
from datetime import datetime,timezone
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from pybhpt.swsh import Yslm
from li_normalized_cloud import LiNormalizedCloud
from li_order4_green import LiOrder4Green
from li_separated_source import LiSeparatedSource,SPINS
from environment_dense_metric import DenseLorenzMetricMode
from environment_metric_sampling import precompute_metric
from environment_source import angular_mode
from environment_cloud import mode_flux
from report_environment_forced_mode import source_grid
from paper_full_tetrad import project_metric
from source_provenance import local_dependency_hashes
ROOT=Path(__file__).resolve().parents[1]

def stamp():return datetime.now(timezone.utc).isoformat()
def enc(z):
    z=np.asarray(z,complex);return np.stack([z.real,z.imag],axis=-1).tolist()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    temp=p.with_suffix('.tmp');temp.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n');temp.replace(p)
def inventory(mc,lmax=5):
    # Equatorially symmetric particle and ell_c=m_c cloud: ell+m even.
    # Static m_g=0 has exactly zero effective orbital flux under real-time freezing.
    return [(ell,m) for ell in range(lmax+1) for m in range(-ell,ell+1)
            if (ell+m)%2==0 and m!=mc and (mc!=2 or m>=-3)]

def weighted_samples(metric,r,theta,a,conjugate=False):
    result=[]
    for t in theta:
        h=metric(float(r),float(t))
        if conjugate:h=h.conjugate()
        v=project_metric(h,r,t,a);g=r+1j*a*np.cos(t);sig=abs(g)**2;d=r*r-2*r+a*a
        result.append([*v[:4],g*v[4],g.conjugate()*v[5],g.conjugate()*v[6],
                       g*v[7],sig*d*v[8],(d*v[8]+v[9])/sig])
    return np.asarray(result)

def spherical_coefficients(metric,r,theta,weights,a,mg,jmax):
    samples=weighted_samples(metric,r,theta,a,conjugate=mg<0)
    Y=np.array([[Yslm(s,j,mg,theta) if j>=max(abs(s),abs(mg)) else np.zeros_like(theta)
                 for s in SPINS] for j in range(jmax+1)])
    return 2*np.pi*np.einsum('jck,kc,k->jc',Y.conjugate(),samples,weights)

def summarize(rows,expected,r0,cloud_ell):
    got={(r['ell'],r['m']) for r in rows}
    if len(got)!=len(rows) or not got.issubset(set(expected)):
        raise ValueError("Repeated/unexpected scalar modes")
    result=dict(completed_modes=len(got),expected_modes=len(expected),
        missing_modes=[list(v) for v in expected if v not in got],
        partial_sums_unit_particle_mass_squared_unit_cloud_mass={
            b:sum(r['flux'][b]['orbital_energy'] for r in rows) for b in ('infinity','horizon')})
    result['complete']=got==set(expected)
    if not result['complete']:return result
    result['total_flux']=result.pop('partial_sums_unit_particle_mass_squared_unit_cloud_mass')
    result['total_flux_Li_units']={b:v/.3**6 for b,v in result['total_flux'].items()}
    reference=ROOT/'docs/environment_reproduction/li_reference_comparison_20260916.json'
    if cloud_ell==1 and reference.exists():
        d=json.loads(reference.read_text())
        result['Li_reference']=dict(file=str(reference.relative_to(ROOT)),sha256=sha(reference),
            type='digitized original vector plot; not author raw arrays')
        result['comparison']=[]
        for ref in d['rows']:
            if ref['rp']!=r0:continue
            b=ref['boundary'];signed=ref['li_magnitude']*(-1 if b=='horizon' else 1)
            value=result['total_flux'][b]
            result['comparison'].append(dict(boundary=b,local=value,Li=signed,
                relative_difference_percent=100*(value/signed-1),
                reference_interpolation_sensitivity=ref['li_log_pchip_vs_linear_relative']))
    return result

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--orbit',type=float,default=20.)
    parser.add_argument('--cloud-ell',type=int,choices=[1,2],default=1)
    parser.add_argument('--L',type=int,default=20);parser.add_argument('--jmax',type=int,default=18)
    parser.add_argument('--metric-q',type=int,default=40)
    parser.add_argument('--nr',type=int,default=12);parser.add_argument('--nh',type=int,default=64)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    if not 1<=args.workers<=4 or args.L<args.jmax+2 or args.metric_q<args.L+8:
        raise ValueError("Need 1..4 workers, spheroidal guard L>=jmax+2 and q>=L+8")
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    paths=['report_li_aligned_flux','li_normalized_cloud','li_order4_green','li_separated_source']
    fingerprints=local_dependency_hashes(ROOT/'src',paths)
    # AST dependency walker does not see the declarative source JSON.
    fingerprints['li_source_terms.json']=sha(ROOT/'src/li_source_terms.json')
    cloud=LiNormalizedCloud(ell=args.cloud_ell)
    grid=SimpleNamespace(orbital_radius=args.orbit,source_inner_offset=5e-4,
        source_outer_radius=320.,radial_order=args.nr,horizon_order=args.nh,horizon_log=True)
    panels,radii,weights=source_grid(grid,cloud)
    x,wt=np.polynomial.legendre.leggauss(args.metric_q);theta=np.arccos(x)
    expected=inventory(cloud.m)
    config=dict(a=.88,alpha=.3,r0=args.orbit,cloud_ell=args.cloud_ell,
        cloud=cloud.provenance,metric_spheroidal_L=args.L,metric_spherical_jmax=args.jmax,
        metric_projection_q=args.metric_q,source_angular_q=64,source_fourier_pmax=12,
        source_fourier_dps=64,source_float_precision='complex128',
        source_formula='covariant-identity corrected Li Appendix',
        source_panels=panels.tolist(),nr=args.nr,nh=args.nh,source_inner_offset=5e-4,
        horizon_green_offset=1e-4,green_infinity_order=4,green_horizon_order=4,
        propagating_green_outer=8000.,bound_green_outer=1000.,radial_rtol=1e-11,
        scalar_lmax=5,scalar_modes=[list(v) for v in expected],
        metric_backend='independent local Lorenz reconstruction, previously checked against protected author modes')
    manifest=out/'execution.json'
    record=dict(status='running',pid=os.getpid(),started_utc=stamp(),config=config,
        implementation_sha256=fingerprints,cloud_mass_audit=cloud.normalization_audit(),completed=[],
        limitations=['Author complex-frequency freeze sequence and BL mass horizon prescription not fully public.',
                    'Metric is independently reconstructed, not author-generated data at these parameters.',
                    'Spherical j<=18 retained from L<=20; guard-mode and angular convergence remain to be tested.',
                    'Finite radial quadrature and source/Green endpoints require convergence.',
                    '64-digit Fourier coefficients do not imply a 64-digit metric/Green pipeline.'])
    if manifest.exists():
        old=json.loads(manifest.read_text())
        if old['config']!=config or old['implementation_sha256']!=fingerprints:
            raise ValueError("Resume configuration/implementation changed")
        record['resumed_from_started_utc']=old['started_utc']
    def persist():
        record['summary']=summarize(record['completed'],expected,args.orbit,args.cloud_ell)
        record['updated_utc']=stamp();save(manifest,record)
    persist()
    try:
        for absmg in sorted({abs(m-cloud.m) for _,m in expected}):
            targets=[v for v in expected if abs(v[1]-cloud.m)==absmg]
            pending=[]
            for ell,m in targets:
                path=out/f'mode_l{ell}_m{m}.json'
                if path.exists():
                    row=json.loads(path.read_text())
                    if row['config']!=config or row['implementation_sha256']!=fingerprints:
                        raise ValueError("Saved mode provenance changed")
                    record['completed'].append(row)
                else:pending.append((ell,m))
            persist()
            if not pending:continue
            record['active_metric_abs_m']=absmg;persist()
            print(f'METRIC |mg|={absmg}, L={args.L}, q={args.metric_q}, {len(radii)} radii',flush=True)
            metric,audit=precompute_metric(DenseLorenzMetricMode(args.orbit,.88,absmg,args.L),
                radii,theta,ROOT/'outputs/metric_cache',args.workers)
            banks={}
            for mg in sorted({m-cloud.m for _,m in pending}):
                banks[mg]=np.array([spherical_coefficients(metric,r,theta,wt,.88,mg,args.jmax) for r in radii])
                bank_path=out/f'metric_spherical_mg{mg}.npz'
                np.savez_compressed(bank_path,radii=radii,coefficients=banks[mg],
                    metadata=json.dumps(dict(config=config,metric_sampling=audit,implementation_sha256=fingerprints),sort_keys=True))
            for ell,m in pending:
                started=time.perf_counter();mg=m-cloud.m
                omega=cloud.omega+mg/(args.orbit**1.5+.88)
                def target(t):return angular_mode(t,ell,m,.88**2*(omega**2-.3**2))[0]
                projector=LiSeparatedSource(a=.88,omega_c=cloud.omega,m_c=cloud.m,
                    metric_m=mg,spherical_lmax=args.jmax,cloud_angular=cloud.angular_state,
                    target_angular=target,quadrature=64,pmax=12,dps=64)
                print(f'SOURCE ell={ell}, m={m}, mg={mg}',flush=True)
                J=np.array([projector.project(r,cloud.radial_state(r)[:,0],coeff)
                            for r,coeff in zip(radii,banks[mg])])
                green=LiOrder4Green(.88,.3,omega,ell,m)
                u=green.insol.sol(radii)[0];v=green.upsol.sol(radii)[0]
                zi=np.sum(weights*u*J)/green.w0;zh=np.sum(weights*v*J)/green.w0
                spread=float(np.max(abs(green.wronskian(radii)/green.w0-1)))
                if not np.all(np.isfinite([zi,zh])) or spread>1e-5:
                    raise ValueError(f'Radial closure failed: W spread={spread}')
                flux=mode_flux(omega,m,cloud.omega,cloud.m,.3,.88,zi,zh)
                cutoffs=[]
                for end in panels[1:]:
                    sel=radii<end
                    ci=np.sum(weights[sel]*u[sel]*J[sel])/green.w0
                    ch=np.sum(weights[sel]*v[sel]*J[sel])/green.w0
                    cutoffs.append(dict(source_outer=float(end),flux=mode_flux(omega,m,cloud.omega,cloud.m,.3,.88,ci,ch)))
                row=dict(status='finite_resolution_mode_completed',ell=ell,m=m,metric_m=mg,
                    omega=omega,config=config,implementation_sha256=fingerprints,
                    metric_sampling=audit,metric_bank_sha256=sha(out/f'metric_spherical_mg{mg}.npz'),
                    source_radii=radii.tolist(),source_weights=weights.tolist(),source=enc(J),
                    horizon_amplitude=enc(zh),up_amplitude=enc(zi),propagating=green.propagating,
                    flux=flux,flux_Li_units={b:{k:v/.3**6 for k,v in values.items()} for b,values in flux.items()},
                    source_outer_cutoff_sequence=cutoffs,wronskian_relative_spread=spread,
                    infinity_boundary=green.boundary_audit,green_outer=green.rmax,
                    completed_utc=stamp(),source_and_green_seconds=time.perf_counter()-started)
                # E=Omega_orbit*L is an independent frequency/convention check.
                om=1/(args.orbit**1.5+.88)
                for b in flux:
                    if abs(flux[b]['orbital_energy']-om*flux[b]['orbital_angular_momentum'])>1e-13*max(abs(flux[b]['orbital_energy']),1e-30):
                        raise ValueError("Orbital energy/angular momentum relation failed")
                save(out/f'mode_l{ell}_m{m}.json',row);record['completed'].append(row);persist()
                print(json.dumps(dict(ell=ell,m=m,flux=flux)),flush=True)
        if not record['summary']['complete']:raise ValueError("Incomplete modal inventory")
        record['status']='complete_finite_resolution_flux_convergence_pending'
        record.pop('active_metric_abs_m',None);record['completed_utc']=stamp();persist()
        print(json.dumps(record['summary'],indent=2),flush=True)
    except Exception as exc:
        record.update(status='failed_resumable',exception=repr(exc),traceback=traceback.format_exc())
        persist();raise
if __name__=='__main__':main()
