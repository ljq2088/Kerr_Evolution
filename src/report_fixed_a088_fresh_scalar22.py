"""Fresh self-consistent a=.88 complex-profile, temporally frozen scalar22.

Independent driver: no production run() background guard is bypassed. Existing
source/grid/metric/Green functions are called explicitly with new provenance.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import time
import traceback
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import brentq, minimize_scalar
from environment_kerr_quasibound_cloud import GeneralKerrQuasiboundCloud
from environment_source import project_source
from environment_lorenz_mode import LorenzMetricMode
from environment_metric_sampling import precompute_metric
from environment_dense_metric import configure_metric_backend
from environment_radial import RadialGreen
from environment_response import SampledResponse
from environment_cloud import mode_flux
from report_environment_forced_mode import source_grid
from source_provenance import source_fingerprint, samples_hash, local_dependency_hashes

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'docs/root_cause_followup_20260917'
BASE=ROOT/'docs/field_alignment_20260917/nonstatic'
FROZEN_MODULE_SHA='3403dbf846cf630f6668fa2ed698272a37c687124a3d193d68b155b43816458b'
def encode(z):return [float(np.real(z)),float(np.imag(z))]
def stamp():return datetime.now(timezone.utc).isoformat()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def tag(r):return f'{r:g}'.replace('.','p')
def save(path,data):
    temp=path.with_suffix('.tmp');temp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');temp.replace(path)
def background_fingerprint(cloud):
    return dict(cloud=cloud.provenance,files=local_dependency_hashes(ROOT/'src',
        ['report_fixed_a088_fresh_scalar22','environment_kerr_quasibound_cloud']))
def finite_bl_integrals(cloud,offset,outer,frequency,nr=12001):
    r=cloud.rp+np.geomspace(offset,outer-cloud.rp,nr);r[0]=cloud.rp+offset;r[-1]=outer
    R,Rp=cloud.radial(r);x,wt=np.polynomial.legendre.leggauss(48)
    S,Sp,_=cloud.angular(np.arccos(x));rr=r[:,None];z=x[None,:]
    sig=rr*rr+cloud.a**2*z*z;d=(rr-cloud.rp)*(rr-cloud.rm)
    gtt=-((rr*rr+cloud.a**2)**2-cloud.a**2*d*(1-z*z))/(sig*d)
    gtp=-2*cloud.a*rr/(sig*d);gpp=(d-cloud.a**2*(1-z*z))/(sig*d*(1-z*z))
    w=complex(frequency);f2=abs(R[:,None])**2*abs(S)**2
    E=sig*((-gtt*abs(w)**2+gpp*cloud.m**2+cloud.mu**2)*f2)
    E+=d*abs(Rp[:,None])**2*abs(S)**2+abs(R[:,None])**2*abs(Sp)**2
    Q=2*sig*(-gtt*w.real+gtp*cloud.m)*f2
    return dict(energy=float(2*np.pi*simpson(E@wt,x=r)),charge=float(2*np.pi*simpson(Q@wt,x=r)))
def normalization_audit(cloud):
    return dict(spectral_KS_T0=cloud.horizon_balance(12001),
        frozen_KS_T0=dict(zip(('energy','charge'),cloud.ks_integrals(12001,frequency=cloud.omega))),
        finite_BL_t0=[dict(inner_offset=e,outer_radius=o,
            spectral=finite_bl_integrals(cloud,e,o,cloud.spectral_omega),
            temporal_frozen=finite_bl_integrals(cloud,e,o,cloud.omega))
            for e,o in [(1e-3,cloud.rmax),(5e-4,cloud.rmax),(1e-5,cloud.rmax),(1e-6,cloud.rmax),(5e-4,320.)]],
        interpretation='Distinct slices and explicit finite BL cutoffs; no mass renormalization is applied. Complex quasibound BL t=0 is not a regular future-horizon-crossing slice. The source-outer320 BL integral also omits the external cloud tail.')
def get_settings(rp):
    path=BASE/f'rp{tag(rp)}_sl2_sm2_L6_q12_nr8_h32.json'
    baseline=json.loads(path.read_text());p=baseline['parameters'];m=p['metric']
    if (m['orbital_radius'],m['m_g'],p['scalar_ell'],p['scalar_m'])!=(rp,1,2,2):raise ValueError('Unexpected baseline')
    rplus=1+np.sqrt(1-m['a']**2)
    inner=p['source_panels'][0]-rplus
    # Baseline stores the subtraction through floating arithmetic. Recover its
    # declared 5e-4 setting only after checking the exact recorded endpoint.
    if abs(inner-5e-4)>1e-14:raise ValueError('Baseline source offset changed')
    args=argparse.Namespace(orbital_radius=rp,source_inner_offset=5e-4,
        source_outer_radius=p['source_panels'][-1],radial_order=p['radial_order'],
        angular_order=p['angular_order'],horizon_order=p['horizon_quadrature_order'],
        horizon_log=p['horizon_log_first_panel'],metric_ellmax=m['ellmax'],
        green_outer_radius=p['green_outer_radius'],green_horizon_offset=p['green_horizon_offset'],
        infinity_method=p['infinity_method'])
    return path,baseline,args

def describe_response(response):
    radii=np.array([3.,10.,20.,40.,50.,80.,100.,150.,200.,320.])
    vals=response.evaluate(radii)
    grid=np.linspace(20,320,1201);f=response.evaluate(grid)[0]
    zero=[];minima=[]
    for i in range(len(grid)-1):
        if f[i].real*f[i+1].real<0:
            zero.append(float(brentq(lambda r:response.evaluate([r])[0,0].real,grid[i],grid[i+1])))
    for i in range(1,len(grid)-1):
        if abs(f[i])<min(abs(f[i-1]),abs(f[i+1])):
            v=minimize_scalar(lambda r:abs(response.evaluate([r])[0,0])**2,
                bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':1e-10})
            minima.append(dict(r=float(v.x),absolute_field=float(np.sqrt(v.fun))))
    return dict(samples=[dict(r=float(r),field=encode(v),derivative=encode(d))
               for r,v,d in zip(radii,*vals)],
        real_part_zeros_no_phase_rotation=zero,absolute_field_local_minima=minima,
        z_h=encode(response.horizon_coefficient),up_coefficient=encode(response.up_coefficient))

def run(rp,workers):
    configure_metric_backend(False)
    baseline_path,baseline,args=get_settings(rp)
    if sha(ROOT/'src/environment_kerr_quasibound_cloud.py')!=FROZEN_MODULE_SHA:
        raise ValueError('New cloud module differs from independently validated immutable version')
    cloud=GeneralKerrQuasiboundCloud(a=.88,alpha=.3,freeze_growth=True)
    provenance=source_fingerprint();background=background_fingerprint(cloud)
    base=LorenzMetricMode(rp,cloud.a,1,args.metric_ellmax)
    panels,radii,weights=source_grid(args,cloud)
    omega=cloud.omega+base.omega
    green=RadialGreen(cloud.a,cloud.mu,omega,2,2,rmax=args.green_outer_radius,
        offset=args.green_horizon_offset,rtol=1e-11,infinity_method=args.infinity_method)
    if not green.rmin<=panels[0]<rp<panels[-1]<=min(cloud.rmax,green.rmax):
        raise ValueError('Source cutoffs outside solved domains')
    p=dict(alpha=.3,cloud_mass=1.,background='general-kerr-quasibound-temporal-frozen',
        metric=base.provenance,scalar_ell=2,scalar_m=2,omega=float(omega),
        radial_order=args.radial_order,angular_order=args.angular_order,
        source_panels=panels.tolist(),green_outer_radius=args.green_outer_radius,
        green_horizon_offset=args.green_horizon_offset,infinity_method=args.infinity_method,
        horizon_log_first_panel=args.horizon_log,horizon_quadrature_order=args.horizon_order)
    path=OUT/f'fresh_a088_rp{tag(rp)}_sl2_sm2_L6_q12_nr8_h32.json'
    samples=[]
    if path.exists():
        old=json.loads(path.read_text())
        for key,expected in [('parameters',p),('source_provenance',provenance),('background_provenance',background)]:
            if old[key]!=expected:raise ValueError(f'Resume {key} differs; refusing to relabel old samples')
        samples=old['samples']
        if old['source_samples_sha256']!=samples_hash(samples):raise ValueError('Resume sample checksum mismatch')
        if len(samples)>len(radii) or any(s['r']!=float(radii[i]) for i,s in enumerate(samples)):
            raise ValueError('Resume source grid differs')
    started=time.perf_counter()
    result=dict(status='source_sampling_in_progress',parameters=p,samples=samples,
        source_provenance=provenance,background_provenance=background,
        baseline_control=dict(file=str(baseline_path.relative_to(ROOT)),sha256=sha(baseline_path),
            settings_read_from_actual_metadata=True,source_cutoffs_same_in_r_minus_rplus=True,
            baseline_a=baseline['parameters']['metric']['a']),
        cloud_normalization_audit=normalization_audit(cloud),
        started_utc=stamp(),pid=os.getpid())
    def persist():
        result['source_samples_sha256']=samples_hash(samples);save(path,result)
    persist()
    theta=np.arccos(np.polynomial.legendre.leggauss(args.angular_order)[0])
    print(f'FRESH A088 ORBIT {rp:g}: metric {len(radii)} radii, workers={workers}',flush=True)
    metric,audit=precompute_metric(base,radii,theta,ROOT/'outputs/metric_cache',workers)
    if source_fingerprint()!=provenance or background_fingerprint(cloud)!=background:
        raise ValueError('Implementation changed during metric sampling')
    result['metric_sampling']=audit;persist()
    for index in range(len(samples),len(radii)):
        frequency,J=project_source(cloud,[float(radii[index])],rp,2,2,metric,ntheta=args.angular_order)
        if abs(frequency-omega)>1e-14:raise ValueError('Source/Green frequency mismatch')
        samples.append(dict(r=float(radii[index]),weight=float(weights[index]),source=encode(J[0]),reused=False))
        persist()
        print(f'FRESH A088 ORBIT {rp:g}: source {index+1}/{len(radii)}',flush=True)
    if source_fingerprint()!=provenance or background_fingerprint(cloud)!=background:
        raise ValueError('Implementation changed during source sampling')
    J=np.array([complex(*s['source']) for s in samples]);u,du=green.insol.sol(radii);v,dv=green.upsol.sol(radii)
    zi=np.sum(weights*u*J)/green.w0;zh=np.sum(weights*v*J)/green.w0
    cutoff=[];fields=[]
    for end in panels[1:]:
        sel=radii<end;ci=np.sum(weights[sel]*u[sel]*J[sel])/green.w0;ch=np.sum(weights[sel]*v[sel]*J[sel])/green.w0
        cutoff.append(dict(outer_source_cutoff=float(end),z_inf=encode(ci if green.propagating else 0j),up_coefficient=encode(ci),z_h=encode(ch)))
    for r in np.concatenate(([green.rmin],panels,[green.rmax])):
        sel=radii<r;cu=np.sum(weights[sel]*u[sel]*J[sel])/green.w0;cv=np.sum(weights[~sel]*v[~sel]*J[~sel])/green.w0
        ur,dur=green.insol.sol(r);vr,dvr=green.upsol.sol(r)
        fields.append(dict(r=float(r),field=encode(vr*cu+ur*cv),derivative=encode(dvr*cu+dur*cv)))
    result.update(status='truncated_single_mode_not_converged',z_inf=encode(zi if green.propagating else 0j),
        up_coefficient=encode(zi),z_h=encode(zh),propagating=bool(green.propagating),
        flux=mode_flux(omega,2,cloud.omega,1,cloud.mu,cloud.a,zi,zh),
        flux_scaling='per q^2*(cloud mass/M); underlying complex cloud normalized to unit KS Killing mass; time growth then frozen',
        radial_response_at_panel_boundaries=fields,outer_source_cutoff_sequence=cutoff,
        wronskian_relative_spread=float(np.max(abs(green.wronskian(radii)/green.w0-1))),
        infinity_boundary_audit=dict(method=green.infinity_method,inverse_r_series_last_term_relative=green.series_last_term_relative,
             status='Coulomb32000 with separate prior convergence controls; no new outer-boundary scan here'),
        completed_utc=stamp(),elapsed_seconds=time.perf_counter()-started,
        missing_convergence=['radial quadrature','horizon source cutoff','outer source cutoff','angular projection','metric ell truncation','other scalar modes'],
        limitations=['Real temporal frequency approximation is explicit; cloud spatial spectrum remains complex.',
             'Fresh metric and cloud share a=.88. No old source is reused or rescaled.',
             'Only scalar22 is recalculated; this does not by itself reproduce full fields or total flux.',
             'The frozen field KG residual is quantified in the separate cloud validation.'])
    persist()
    fresh_response=SampledResponse.from_report(result);old_response=SampledResponse.from_report(baseline)
    new=describe_response(fresh_response);old=describe_response(old_response)
    comparison=dict(status='fresh_a088_versus_threshold_same_numerical_settings',
        fresh_file=str(path.relative_to(ROOT)),fresh_sha256=sha(path),baseline_file=str(baseline_path.relative_to(ROOT)),baseline_sha256=sha(baseline_path),
        fresh=new,baseline=old,quadrature_flux_fresh=result['flux'],quadrature_flux_baseline=baseline['flux'],
        radial_comparison=[dict(r=n['r'],new_field=n['field'],old_field=o['field'],
            amplitude_ratio=abs(complex(*n['field']))/abs(complex(*o['field'])),
            complex_relative_change=encode((complex(*n['field'])-complex(*o['field']))/complex(*o['field'])))
            for n,o in zip(new['samples'],old['samples'])],
        z_h_complex_relative_change=encode((fresh_response.horizon_coefficient-old_response.horizon_coefficient)/old_response.horizon_coefficient),
        nominal_up_coefficient_comparison_caution='Bound Up has arbitrary endpoint normalization and depends strongly on kappa; compare physical field instead of coefficient magnitudes.',
        no_fitted_amplitude_or_phase=True,
        limitation='This changes physical background a, cloud complex eigenfunction and KS normalization consistently, at a fixed finite numerical resolution; it is not a convergence bound or full mode-summed field.')
    cpath=OUT/f'fresh_a088_rp{tag(rp)}_scalar22_comparison.json';save(cpath,comparison)
    print(json.dumps(dict(orbit=rp,result=str(path),comparison=str(cpath),flux=result['flux']),indent=2),flush=True)
    return dict(orbital_radius=rp,file=str(path.relative_to(ROOT)),sha256=sha(path),comparison=str(cpath.relative_to(ROOT)),elapsed_seconds=result['elapsed_seconds'])

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--orbits',type=float,nargs='+',default=[42.1,41.1]);args=parser.parse_args()
    if args.workers<1 or args.workers>4:raise ValueError('Bounded run permits 1..4 radius workers')
    OUT.mkdir(parents=True,exist_ok=True)
    manifest=OUT/'fresh_a088_scalar22_execution.json'
    record=dict(status='running',started_utc=stamp(),pid=os.getpid(),workers=args.workers,
        orbits=args.orbits,driver_sha256=sha(__file__),cloud_module_sha256=FROZEN_MODULE_SHA,results=[])
    save(manifest,record)
    try:
        for rp in args.orbits:
            record['active_orbit']=rp;save(manifest,record)
            record['results'].append(run(rp,args.workers));save(manifest,record)
        record['status']='finite_resolution_fresh_scalar22_completed';record.pop('active_orbit',None)
        record['completed_utc']=stamp();save(manifest,record)
    except Exception as exc:
        record.update(status='failed_resumable',exception=repr(exc),traceback=traceback.format_exc(),failed_utc=stamp())
        save(manifest,record);raise
if __name__=='__main__':main()
