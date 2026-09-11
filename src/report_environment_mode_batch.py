"""Solve several scalar multipoles using one cached nonstatic metric mode.

Every scalar channel retains its own source projection, radial operator,
response file, and convergence requirements. Only the metric/cloud are shared.
"""
import argparse
import hashlib
import json
from pathlib import Path
from time import perf_counter
from report_environment_forced_mode import argument_parser,build_cloud,run
from environment_lorenz_mode import LorenzMetricMode,ConjugateMetricMode


def main():
    parser=argparse.ArgumentParser(add_help=False)
    parser.add_argument('--scalar-ells',nargs='+',type=int,required=True)
    parser.add_argument('--conjugate-ells',nargs='+',type=int,default=[])
    batch,remaining=parser.parse_known_args()
    args=argument_parser().parse_args(remaining)
    if args.scalar_ell is not None or args.reuse_source or args.extend_metric_source:
        raise ValueError('Specify --scalar-ells; scalar-specific reuse/extension belongs in the single-channel command')
    if len(set(batch.scalar_ells))!=len(batch.scalar_ells) or len(set(batch.conjugate_ells))!=len(batch.conjugate_ells):
        raise ValueError('Duplicate scalar multipoles')
    cloud=build_cloud(args)
    metric=LorenzMetricMode(args.orbital_radius,cloud.a,args.metric_m,args.metric_ellmax)
    if min(batch.scalar_ells)<abs(args.metric_m+cloud.m):
        raise ValueError('Every scalar ell must be >= |m_g+m_cloud|')
    if batch.conjugate_ells and min(batch.conjugate_ells)<abs(cloud.m-args.metric_m):
        raise ValueError('Opposite metric branch requires ell >= |m_cloud-m_g|')
    parameters=dict(vars(args),scalar_ells=batch.scalar_ells,conjugate_ells=batch.conjugate_ells)
    fingerprint=hashlib.sha256(json.dumps(parameters,sort_keys=True,default=str).encode()).hexdigest()[:12]
    folder=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    output=folder/f'scalar_batch_mg{args.metric_m}_{fingerprint}.json'
    summary=dict(status='batch_in_progress',parameters=parameters,channels=[])
    def save():
        temporary=output.with_suffix('.tmp')
        temporary.write_text(json.dumps(summary,indent=2,default=str)+'\n')
        temporary.replace(output)
    save()
    plan=[(ell,metric,False) for ell in batch.scalar_ells]
    if batch.conjugate_ells:
        conjugate=ConjugateMetricMode(metric)
        plan.extend((ell,conjugate,True) for ell in batch.conjugate_ells)
    for ell,selected_metric,conjugated in plan:
        selected=argparse.Namespace(**vars(args));selected.scalar_ell=ell;selected.metric_m=selected_metric.m
        before=metric._values.cache_info();started=perf_counter()
        path,result=run(selected,cloud=cloud,metric=selected_metric)
        after=metric._values.cache_info()
        channel=dict(scalar_ell=ell,scalar_m=cloud.m+selected_metric.m,file=path.name,
            metric_m=selected_metric.m,metric_conjugated=conjugated,
            metric_cache_hits=after.hits-before.hits,metric_reconstructions=after.misses-before.misses,
            elapsed_seconds=perf_counter()-started,flux=result['flux'])
        summary['channels'].append(channel);save()
        print('Completed shared-metric channel',json.dumps(channel),flush=True)
    summary['status']='batch_completed_finite_resolution_not_converged';save()
    print(output.name,flush=True)


if __name__=='__main__':main()
