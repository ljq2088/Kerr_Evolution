"""Solve several scalar multipoles using one cached nonstatic metric mode.

Every scalar channel retains its own source projection, radial operator,
response file, and convergence requirements. Only the metric/cloud are shared.
"""
import argparse
import hashlib
import json
from pathlib import Path
from time import perf_counter
from report_environment_forced_mode import argument_parser,build_cloud,build_metric,source_grid,run
from environment_lorenz_mode import LorenzMetricMode,ConjugateMetricMode


def main():
    parser=argparse.ArgumentParser(add_help=False)
    parser.add_argument('--scalar-ells',nargs='+',type=int,default=[])
    parser.add_argument('--conjugate-ells',nargs='+',type=int,default=[])
    parser.add_argument('--workers',type=int,default=1)
    batch,remaining=parser.parse_known_args()
    args=argument_parser().parse_args(remaining)
    if args.scalar_ell is not None or args.reuse_source or args.extend_metric_source or args.output:
        raise ValueError('Specify --scalar-ells; scalar-specific reuse/extension belongs in the single-channel command')
    if len(set(batch.scalar_ells))!=len(batch.scalar_ells) or len(set(batch.conjugate_ells))!=len(batch.conjugate_ells):
        raise ValueError('Duplicate scalar multipoles')
    if not batch.scalar_ells and not batch.conjugate_ells:
        raise ValueError('At least one direct or conjugate scalar multipole is required')
    if batch.workers<1:raise ValueError('Workers must be positive')
    cloud=build_cloud(args)
    metric=build_metric(args,cloud)
    if metric.m==0 and batch.conjugate_ells:
        raise ValueError('Static metric has no distinct opposite-m branch')
    if batch.scalar_ells and min(batch.scalar_ells)<abs(args.metric_m+cloud.m):
        raise ValueError('Every scalar ell must be >= |m_g+m_cloud|')
    if batch.conjugate_ells and min(batch.conjugate_ells)<abs(cloud.m-args.metric_m):
        raise ValueError('Opposite metric branch requires ell >= |m_cloud-m_g|')
    parameters=dict(vars(args),scalar_ells=batch.scalar_ells,conjugate_ells=batch.conjugate_ells)
    parameters.pop('output',None)
    if not parameters.get('dense_angular'):
        parameters.pop('dense_angular',None)
    fingerprint=hashlib.sha256(json.dumps(parameters,sort_keys=True,default=str).encode()).hexdigest()[:12]
    folder=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    output=folder/f'scalar_batch_mg{args.metric_m}_{fingerprint}.json'
    summary=dict(status='batch_in_progress',parameters=parameters,channels=[])
    def save():
        temporary=output.with_suffix('.tmp')
        temporary.write_text(json.dumps(summary,indent=2,default=str)+'\n')
        temporary.replace(output)
    save()
    if batch.workers<1:raise ValueError('Workers must be positive')
    if batch.workers>1:
        import numpy as np
        from environment_metric_sampling import precompute_metric
        _,radii,_=source_grid(args,cloud)
        theta=np.arccos(np.polynomial.legendre.leggauss(args.angular_order)[0])
        try:
            metric,audit=precompute_metric(metric,radii,theta,folder.parents[1]/'outputs/metric_cache',batch.workers)
        except Exception as error:
            summary.update(status='batch_failed_partial_results_preserved',
                failure=dict(stage='metric_precomputation',error_type=type(error).__name__,message=str(error)))
            save()
            raise
        summary['metric_precomputation']=audit;save()
    plan=[(ell,metric,False) for ell in batch.scalar_ells]
    if batch.conjugate_ells:
        conjugate=ConjugateMetricMode(metric)
        plan.extend((ell,conjugate,True) for ell in batch.conjugate_ells)
    for ell,selected_metric,conjugated in plan:
        selected=argparse.Namespace(**vars(args));selected.scalar_ell=ell;selected.metric_m=selected_metric.m
        before=metric._values.cache_info();started=perf_counter()
        try:
            path,result=run(selected,cloud=cloud,metric=selected_metric)
        except Exception as error:
            summary.update(status='batch_failed_partial_results_preserved',
                failure=dict(scalar_ell=ell,scalar_m=cloud.m+selected_metric.m,
                             error_type=type(error).__name__,message=str(error)))
            save()
            raise
        after=metric._values.cache_info()
        channel=dict(scalar_ell=ell,scalar_m=cloud.m+selected_metric.m,file=path.name,
            metric_m=selected_metric.m,metric_conjugated=conjugated,
            metric_cache_hits=after.hits-before.hits+(after.misses-before.misses if getattr(metric,'precomputed',False) else 0),
            metric_reconstructions=0 if getattr(metric,'precomputed',False) else after.misses-before.misses,
            elapsed_seconds=perf_counter()-started,flux=result['flux'])
        summary['channels'].append(channel);save()
        print('Completed shared-metric channel',json.dumps(channel),flush=True)
    summary['status']='batch_completed_finite_resolution_not_converged';save()
    print(output.name,flush=True)


if __name__=='__main__':main()
