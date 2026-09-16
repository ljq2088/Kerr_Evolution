"""Parallel fresh source calculation for the remaining H00 convergence audit."""
import argparse,json
from pathlib import Path
import numpy as np
from report_environment_forced_mode import argument_parser,build_cloud,build_metric,source_grid,run
from environment_metric_sampling import precompute_metric
from source_provenance import source_fingerprint,validate_saved_samples


def main():
    p=argparse.ArgumentParser();p.add_argument('--radial-order',type=int,required=True);p.add_argument('--horizon-order',type=int,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--workers',type=int,default=4);choice=p.parse_args()
    args=argument_parser().parse_args(['--output',str(choice.output),'--radial-order',str(choice.radial_order),'--angular-order','10','--metric-ellmax','4','--metric-m','-1','--scalar-ell','0','--source-inner-offset','.0005','--source-outer-radius','320','--horizon-log','--horizon-order',str(choice.horizon_order)])
    cloud=build_cloud(args);metric=build_metric(args,cloud)
    _,radii,_=source_grid(args,cloud);theta=np.arccos(np.polynomial.legendre.leggauss(10)[0]);completed=0
    if choice.output.exists():
        old=json.loads(choice.output.read_text());validate_saved_samples(old,source_fingerprint());completed=len(old['samples'])
    metric,audit=precompute_metric(metric,radii[completed:],theta,Path(__file__).resolve().parents[1]/'outputs/metric_cache',choice.workers)
    out,result=run(args,cloud=cloud,metric=metric)
    result['execution_audit']=dict(metric_precomputation=audit,already_completed_same_provenance_sources=completed,note='Sequential sampler was interrupted only to use bounded parallel metric precomputation. All existing source samples had the same recorded source fingerprint; no unversioned historical source was imported.')
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
