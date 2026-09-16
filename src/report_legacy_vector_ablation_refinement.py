"""Vector-only refinement; do not mix refined vector and coarse full integrals."""
import argparse,json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
from types import SimpleNamespace
import numpy as np
from report_legacy_vector_ablation import VectorDipole,enc,sha,ROOT
from environment_source import ThresholdCloud,project_source
from environment_radial import RadialGreen
from report_environment_forced_mode import source_grid

CLOUD=None
METRIC=None
ORBIT=None
NQ=None
def one_radius(r):
    _,J=project_source(CLOUD,[r],ORBIT,0,0,METRIC,ntheta=NQ)
    return J[0]

def main():
    global CLOUD,METRIC,ORBIT,NQ
    p=argparse.ArgumentParser();p.add_argument('--coarse',required=True);p.add_argument('--output',required=True);p.add_argument('--radial-order',type=int,default=16);p.add_argument('--horizon-order',type=int,default=64);p.add_argument('--angular-order',type=int,default=24);p.add_argument('--jet-order',type=int,default=8);p.add_argument('--workers',type=int,default=1);args=p.parse_args()
    previous=json.loads(Path(args.coarse).read_text());
    if previous['status']!='completed_same_grid_vector_ablation':raise ValueError('Coarse experiment must be complete')
    par=previous['parameters'];cloud=ThresholdCloud(alpha=par['alpha']);orbit=par['metric']['orbital_radius']
    gridargs=SimpleNamespace(source_inner_offset=par['source_panels'][0]-cloud.rp,source_outer_radius=par['source_panels'][-1],orbital_radius=orbit,radial_order=args.radial_order,horizon_order=args.horizon_order,horizon_log=par['horizon_log_first_panel'])
    panels,rs,weights=source_grid(gridargs,cloud);omega=cloud.omega-1/(orbit**1.5+cloud.a)
    green=RadialGreen(cloud.a,cloud.mu,omega,0,0,rmax=par['green_outer_radius'],offset=par['green_horizon_offset'],rtol=5e-13)
    kernel=green.upsol.sol(rs)[0]/green.w0;metric=VectorDipole(orbit,cloud.a,args.jet_order)
    out=Path(args.output)
    if out.exists():raise FileExistsError(out)
    result=dict(status='sampling_refined_vector_only',coarse=str(args.coarse),coarse_sha256=sha(args.coarse),parameters=par,refined_controls=dict(radial_order=args.radial_order,horizon_order=args.horizon_order,angular_order=args.angular_order,jet_order=args.jet_order,green_rtol=5e-13,workers=args.workers),scope='Only Z_vector refinement. No refined full source exists on this new grid, so no mixed-grid legacy or physical flux is reported.',samples=[])
    start=time.perf_counter()
    def save():
        temp=out.with_suffix('.tmp');temp.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');temp.replace(out)
    save()
    CLOUD,METRIC,ORBIT,NQ=cloud,metric,orbit,args.angular_order
    with ProcessPoolExecutor(max_workers=args.workers,mp_context=mp.get_context('fork')) as pool:
        for i,Jvalue in enumerate(pool.map(one_radius,rs,chunksize=1)):
            r=rs[i]
            result['samples'].append(dict(r=float(r),weight=float(weights[i]),source_vector=enc(Jvalue),green_up_over_W=enc(kernel[i])))
            if (i+1)%16==0 or i+1==len(rs):save();print(f'refined vector {i+1}/{len(rs)}',flush=True)
    J=np.array([complex(*x['source_vector']) for x in result['samples']]);Z=np.dot(weights,kernel*J);old=complex(*previous['z_h_vector'])
    result.update(status='completed_refined_vector_only',z_h_vector=enc(Z),coarse_z_h_vector=enc(old),relative_Z_vector_change=float(abs(Z-old)/abs(Z)),absolute_Z_vector_change=float(abs(Z-old)),elapsed_seconds=time.perf_counter()-start)
    save();print('Zvector refined',Z,'coarse',old,'relative',result['relative_Z_vector_change'],flush=True)
if __name__=='__main__':main()
