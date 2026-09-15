"""Fixed-source precision audit of all 21 alpha=.3 rp=20 massive responses.

Only homogeneous scalar radial algorithms change: historical source samples and
finite quadrature are held fixed, with explicit input hashes. This is not a
reconstruction or source convergence audit. Massive KG is not the massless
Teukolsky equation, so MST/GSN cannot be substituted unchanged.
"""
import json, hashlib, time
from contextlib import contextmanager
from pathlib import Path
import numpy as np
import environment_radial as radial
from environment_cloud import mode_flux

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'docs/environment_reproduction/massive_radial_precision_audit.json'

@contextmanager
def boundary_order(order):
    original=radial.infinity_series
    def expanded(*args, **kwargs):
        # RadialGreen requests order six and five for its error indicator.
        requested=kwargs.pop('order',6)
        return original(*args,order=order if requested==6 else order-1,**kwargs)
    radial.infinity_series=expanded
    try:
        yield
    finally:
        radial.infinity_series=original

def encode(z): return [float(z.real),float(z.imag)]

def calculate(data,config):
    p=data['parameters']; a=p['metric']['a']; mu=p['alpha']
    ell,m,w=p['scalar_ell'],p['scalar_m'],p['omega']
    radii=np.array([s['r'] for s in data['samples']])
    weights=np.array([s['weight'] for s in data['samples']])
    source=np.array([complex(*s['source']) for s in data['samples']])
    with boundary_order(config['order']):
        g=radial.RadialGreen(a,mu,w,ell,m,rmax=config['rmax'],offset=config['offset'],
                            rtol=config['rtol'],infinity_method=config['method'])
    u,du=g.insol.sol(radii); v,dv=g.upsol.sol(radii)
    kh,ki=v/g.w0,u/g.w0
    zh=np.sum(weights*source*kh); zi=np.sum(weights*source*ki) if g.propagating else 0j
    omega_c=a/(2*g.rp)
    flux=mode_flux(w,m,omega_c,1,mu,a,zi,zh)
    aa=u*dv; bb=v*du
    condition=(abs(aa)+abs(bb))/np.maximum(abs(aa-bb),1e-300)
    summary=dict(z_h=encode(zh),z_inf=encode(zi),
                 horizon_flux=flux['horizon']['orbital_energy'],
                 infinity_flux=flux['infinity']['orbital_energy'],
                 wronskian_spread=float(np.max(abs(g.wronskian(radii)/g.w0-1))),
                 wronskian_max_subtraction_condition=float(max(condition)),
                 horizon_integral_cancellation=float(np.sum(abs(weights*source*kh))/max(abs(zh),1e-300)),
                 infinity_integral_cancellation=float(np.sum(abs(weights*source*ki))/max(abs(zi),1e-300)) if g.propagating else None,
                 boundary_last_term_relative=g.series_last_term_relative,
                 k=encode(complex(g.k)),horizon_frequency=float(w-m*g.oh),
                 beta_over_k2_rmax=float(abs((2*w*w-mu*mu)/(g.k*g.k*config['rmax']))),
                 max_In_magnitude=float(max(abs(u))),max_Up_magnitude=float(max(abs(v))),
                 In_evaluations=g.insol.nfev,Up_evaluations=g.upsol.nfev)
    return summary,kh,ki

def main():
    audit=json.loads((ROOT/'docs/environment_reproduction/full_reference_audit_20260915.json').read_text())
    names=[n for n in audit['input_response_sha256'] if n.startswith('forced_mode_nr8_nt18_L18_')
           and 'alpha' not in n and 'schwarzschild' not in n]
    configs={
      'baseline':dict(rmax=1000.,offset=1e-4,rtol=1e-11,order=6,method='series'),
      'tight_tolerance':dict(rmax=1000.,offset=1e-4,rtol=5e-14,order=6,method='series'),
      'boundary_refined':dict(rmax=2000.,offset=1e-5,rtol=1e-12,order=12,method='series'),
      'coulomb_propagating':dict(rmax=4000.,offset=1e-5,rtol=1e-12,order=12,method='coulomb'),
    }
    results=dict(status='running_fixed_source_algorithm_precision_diagnostic',configurations=configs,rows=[],
       limitations=['Source arrays and finite quadrature remain shared; no claim of independent metric/source validation.',
                    'Coulomb has only r^-2 far potential; rmax convergence remains necessary.',
                    'Relative changes are dimensionless ratios, not percent.'],
       implementation_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()
         for f in [Path(__file__),ROOT/'src/environment_radial.py',ROOT/'src/environment_cloud.py']})
    def save():
        tmp=OUT.with_suffix('.tmp');tmp.write_text(json.dumps(results,indent=2)+'\n');tmp.replace(OUT)
    for filename in names:
        path=OUT.parent/filename;data=json.loads(path.read_text());p=data['parameters']
        row=dict(file=filename,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
          ell=p['scalar_ell'],m=p['scalar_m'],omega=p['omega'],metrics={})
        base_h=base_i=None
        for name,config in configs.items():
            if name=='coulomb_propagating' and not data['propagating']:continue
            t=time.perf_counter();s,kh,ki=calculate(data,config)
            if name=='baseline':
                base_h,base_i=kh,ki
                s['stored_horizon_flux_relative_change']=s['horizon_flux']/data['flux']['horizon']['orbital_energy']-1
            else:
                s['horizon_flux_relative_change']=s['horizon_flux']/row['metrics']['baseline']['horizon_flux']-1
                s['infinity_flux_relative_change']=(s['infinity_flux']/row['metrics']['baseline']['infinity_flux']-1 if data['propagating'] else None)
                s['kh_max_relative_change']=float(max(abs(kh/base_h-1)))
                if data['propagating']:s['ki_max_relative_change']=float(max(abs(ki/base_i-1)))
            s['seconds']=time.perf_counter()-t
            row['metrics'][name]=s
            print(json.dumps(dict(ell=row['ell'],m=row['m'],config=name,seconds=s['seconds'],
                 fh=s['horizon_flux'],fi=s['infinity_flux'],Wspread=s['wronskian_spread'],
                 dH=s.get('horizon_flux_relative_change'),dI=s.get('infinity_flux_relative_change'))),flush=True)
        results['rows'].append(row);save()
    totals={}
    for name in configs:
        if name=='coulomb_propagating':continue
        fh=sum(r['metrics'][name]['horizon_flux'] for r in results['rows'])
        fi=sum(r['metrics'][name]['infinity_flux'] for r in results['rows'])
        totals[name]=dict(horizon_flux=fh,infinity_flux=fi)
    for name in totals:
        if name!='baseline':
            for boundary in ['horizon','infinity']:
                totals[name][boundary+'_relative_change']=totals[name][boundary+'_flux']/totals['baseline'][boundary+'_flux']-1
    # Coulomb applies only to propagating channels; replace those nine contributions.
    for boundary in ['horizon','infinity']:
        total=sum(r['metrics'].get('coulomb_propagating',r['metrics']['baseline'])[boundary+'_flux'] for r in results['rows'])
        totals.setdefault('coulomb_propagating_mixed',{})[boundary+'_flux']=total
        totals['coulomb_propagating_mixed'][boundary+'_relative_change']=total/totals['baseline'][boundary+'_flux']-1
    coverage_path=OUT.parent/'flux_coverage_L18_nt18_i6_h5_f5.json'
    coverage=json.loads(coverage_path.read_text())
    byfile={r['file']:r for r in results['rows']}
    coverage_totals={}
    for name in configs:
        block={}
        for boundary in ['horizon','infinity']:
            block[boundary+'_flux']=sum(byfile[c['file']]['metrics'].get(name,byfile[c['file']]['metrics']['baseline'])[boundary+'_flux']
                                       for c in coverage[boundary]['modes'])
        coverage_totals[name]=block
    for name,block in coverage_totals.items():
        if name!='baseline':
            for boundary in ['horizon','infinity']:
                block[boundary+'_relative_change']=block[boundary+'_flux']/coverage_totals['baseline'][boundary+'_flux']-1
    results.update(status='completed_fixed_source_algorithm_precision_diagnostic',
       totals_all_21_both_boundaries=totals,paper_coverage_totals=coverage_totals,
       coverage_sha256=hashlib.sha256(coverage_path.read_bytes()).hexdigest(),
       coverage_note='Original published-comparison report includes only ell<=5 at horizon, ell<=6 at infinity; all-21 sum additionally includes scalar (6,2) horizon flux.')
    save();print(json.dumps(coverage_totals,indent=2))

if __name__=='__main__':main()
