"""Recompute actual scalar sources/responses with versioned output files.

Use only metric tensors accepted by the existing code-hashed cache; never
reuse source samples. Compare fresh responses with historical finite reports.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from report_environment_forced_mode import argument_parser,build_cloud,build_metric,source_grid,run
from environment_metric_sampling import precompute_metric
from environment_lorenz_mode import ConjugateMetricMode


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--metric-m',type=int,default=1)
    parser.add_argument('--ellmax',type=int,default=18);parser.add_argument('--workers',type=int,default=4)
    choice=parser.parse_args();mg=choice.metric_m
    baseargs=['--metric-m',str(mg),'--metric-ellmax',str(choice.ellmax),'--radial-order','8',
       '--angular-order','18','--source-inner-offset','.0005','--source-outer-radius','320',
       '--horizon-log','--horizon-order','32']
    args=argument_parser().parse_args(baseargs)
    cloud=build_cloud(args);metric=build_metric(args,cloud)
    _,radii,_=source_grid(args,cloud);theta=np.arccos(np.polynomial.legendre.leggauss(18)[0])
    root=Path(__file__).resolve().parents[1];folder=root/'docs/environment_reproduction'
    # Inspect cache coverage before the explicit precomputation step.
    metric,audit=precompute_metric(metric,radii,theta,root/'outputs/metric_cache',choice.workers)
    rows=[]
    choices=[(ell,metric) for ell in range(abs(1+mg),7) if (ell+1+mg)%2==0]
    opposite=ConjugateMetricMode(metric)
    choices += [(ell,opposite) for ell in range(abs(1-mg),6) if (ell+1-mg)%2==0]
    if choice.ellmax==1:choices=[(0,opposite)]
    for ell,selected in choices:
        current=argparse.Namespace(**vars(args));current.scalar_ell=ell;current.metric_m=selected.m
        current.output=folder/f'fresh_20260915_L{choice.ellmax}_mg{selected.m}_sl{ell}.json'
        path,result=run(current,cloud=cloud,metric=selected)
        suffix='' if selected.m==2 and ell==3 else f'_mg{selected.m}_sl{ell}'
        oldpath=folder/f'forced_mode_nr8_nt18_L{choice.ellmax}{suffix}_inner0.0005_outer320_log_h32.json'
        row=dict(scalar_ell=ell,scalar_m=1+selected.m,new_file=path.name,metric_cache=audit)
        if oldpath.exists():
            old=json.loads(oldpath.read_text());oldJ=np.array([complex(*s['source']) for s in old['samples']]);newJ=np.array([complex(*s['source']) for s in result['samples']])
            if result['parameters']!=old['parameters']:raise ValueError('Physics mismatch')
            np.testing.assert_array_equal([s['r'] for s in result['samples']],[s['r'] for s in old['samples']])
            row.update(old_file=oldpath.name,old_sha256=hashlib.sha256(oldpath.read_bytes()).hexdigest(),
              relative_source_max=float(np.max(abs(newJ-oldJ))/np.max(abs(oldJ))),
              horizon_flux_relative_change=float(result['flux']['horizon']['orbital_energy']/old['flux']['horizon']['orbital_energy']-1),
              infinity_flux_relative_change=(float(result['flux']['infinity']['orbital_energy']/old['flux']['infinity']['orbital_energy']-1) if old['flux']['infinity']['orbital_energy'] else None))
        rows.append(row)
        summary=dict(status='fresh_source_response_comparison_not_paper_reproduction',rows=rows)
        output=folder/f'fresh_environment_audit_L{choice.ellmax}_mg{mg}.json'
        output.write_text(json.dumps(summary,indent=2)+'\n')
        print('COMPARISON',json.dumps(row),flush=True)


if __name__=='__main__':main()
