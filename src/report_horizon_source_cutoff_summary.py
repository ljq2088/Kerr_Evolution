"""Summarize fresh horizon cutoff diagnostics without regenerating physics."""
import hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'docs/environment_reproduction'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def z(v):return complex(*v)
def enc(v):return [float(v.real),float(v.imag)]
def main():
    names=['horizon_source_cutoff_r20_domain_extended.json','horizon_source_cutoff_r20_extended_precision.json','horizon_source_cutoff_r20_failed_controls.json','fresh_20260915_L18_mg-1_sl0.json','metric_backend_response_L1_s2AUTO_gAUTO_rtolNone.json']
    data={n:json.loads((D/n).read_text()) for n in names}
    standard,fine=[data[n] for n in names[:2]]
    baseline=data[names[3]];dipole=data[names[4]]
    from source_provenance import source_fingerprint,validate_saved_samples
    validate_saved_samples(baseline,source_fingerprint())
    for module,digest in dipole['implementation_sha256'].items():
        path=ROOT/'src'/f'{module}.py'
        if path.is_file() and sha(path)!=digest:raise ValueError(f'Stale dipole dependency: {module}')
    rows=[]
    for s,f in zip(standard['rows'],fine['rows']):
        assert s['offset']==f['offset']
        rows.append(dict(offset=f['offset'],delta_Z=f['delta_Z'],relative_flux_change=f['relative_flux_change'],precision_check_delta_Z_absolute=abs(z(f['delta_Z'])-z(s['delta_Z'])),precision_check_relative_flux_change_difference=f['relative_flux_change']-s['relative_flux_change']))
    fits=[f for f in fine['fits'] if f['sample_count']>1+2*f['degree']]
    fitrange=[min(f['relative_flux_change'] for f in fits),max(f['relative_flux_change'] for f in fits)]
    near=[]
    for f,l in zip(baseline['samples'][:8],dipole['samples'][:8]):
        assert f['r']==l['r'] and f['weight']==l['weight']
        delta=z(f['source'])-z(l['source'])
        near.append(dict(r=f['r'],full_L18_source=f['source'],dipole_source=l['source'],source_difference=enc(delta),difference_relative_to_full=abs(delta)/abs(z(f['source']))))
    result=dict(status='completed_bounded_horizon_source_audit',parameters= fine['parameters'],baseline=dict(file=names[3],z_h=baseline['z_h'],flux=baseline['flux']['horizon']['orbital_energy'],provenance_validated=True),references={n:dict(sha256=sha(D/n)) for n in names},implementation=dict(report_horizon_source_cutoff=sha(ROOT/'src/report_horizon_source_cutoff.py'),report_horizon_source_cutoff_summary=sha(Path(__file__))),comparison=rows,accepted_overdetermined_fits=fits,estimated_limit_relative_flux_change_range=fitrange,fit_error_interpretation='Fit spread/residual and two precision configurations are empirical consistency indicators, not a rigorous error bound. The degree2 five-point fit has five coefficients, so its exact interpolation residual is explicitly excluded.',near_existing_cutoff_high_ell_context=dict(angular_note='Full source uses original pybhpt projection; backend L1 audit used independently checked DenseRealHarmonic. Same actual nodes and weights; all recorded local code hashes match.',samples=near,interpretation='Observed differences at radii above the old cutoff are small, but they do not bound high-ell source errors below that cutoff or prove the complete L18 horizon tail converged.'),failed_control=data[names[2]],domain_hazard=dict(function='lorenz_kappa.trace_field_jet',line=33,origin='scalar_resolvent creates RadialGreen(offset=1e-4) at line21; scipy OdeSolution silently extrapolates when called at smaller r.',actual_effect='Uncontrolled future cutoff tests can sample trace/kappa outside solved domains. Existing baseline cutoff5e-4 remains inside the domain.',diagnostic_mitigation='Only process-local lorenz_kappa.RadialGreen constructor binding is replaced; cached resolvents are cleared; integration begins at r+1e-7 and all evaluation radii exceed cloud start r+1e-6.',minimal_fix='Before evaluating radial.insol.sol or radial.upsol.sol in trace_field_jet, reject nonfinite r or r outside [radial.rmin,radial.rmax]. Both trace and finite-difference kappa use this function.',regression_cases=['Mock a resolver with known [rmin,rmax], assert values just below and above raise ValueError before its dense-output callable is invoked.','Assert nan and infinities are rejected.','Both exact endpoints and representative interior points must return the unchanged valid trace jet; verify both r<r0 and r>=r0 branches.','Through kappa_jet confirm invalid g.r raises before any out-of-domain resolvent evaluation.','Separately run a valid current actual alpha.3 rp20 trace/kappa point and compare complex values before/after guard.'],optional_followup='Expose explicit horizon offset through scalar_resolvent, trace_field_jet and kappa_jet with the offset in the cache key; avoid silently changing or automatically extending the solver.'),conclusions=['The dominant scalar00 dipole tail omitted by cutoff5e-4 changes the current fullL18 baseline flux by about -0.0454%, not several percent.','This includes oscillatory Frobenius branches through the complex Green endpoint identity; no orbital term or arbitrary prefactor was added.','No production numerical file was modified.','This is a dipole shell correction on a fullL18 baseline, not an all-ell horizon source convergence certificate.','The newly reported literature horizon normalization correction is separate; no such correction is applied by this diagnostic.'])
    out=D/'horizon_source_cutoff_r20_audit.json';out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(out)
if __name__=='__main__':main()
