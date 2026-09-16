"""Add bounded parallel workers without changing a running metric sampler.

The auxiliary sampler owns a separate cache directory. Only fully validated
atomic cache records with identical metadata are copied to the standard cache.
No process is stopped by this tool.
"""
import argparse,hashlib,json,os,shutil
from pathlib import Path
import numpy as np
from report_environment_forced_mode import argument_parser,build_cloud,build_metric,source_grid
from environment_metric_sampling import precompute_metric
ROOT=Path(__file__).resolve().parents[1]

def main():
 parser=argparse.ArgumentParser(add_help=False);parser.add_argument('--workers',type=int,default=4);parser.add_argument('--output',type=Path,required=True);extra,remaining=parser.parse_known_args();args=argument_parser().parse_args(remaining)
 cloud=build_cloud(args);base=build_metric(args,cloud);_,radii,_=source_grid(args,cloud);theta=np.arccos(np.polynomial.legendre.leggauss(args.angular_order)[0]);selected=radii[len(radii)//2:]
 record=dict(status='auxiliary_tail_sampling',radius=args.orbital_radius,worker_count=extra.workers,all_radii=len(radii),tail_radii=len(selected),driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 extra.output.write_text(json.dumps(record,indent=2)+'\n')
 _,audit=precompute_metric(base,selected,theta,ROOT/'outputs/metric_cache_rp15_auxiliary',extra.workers)
 source=ROOT/'outputs/metric_cache_rp15_auxiliary'/audit['cache_key'];target=ROOT/'outputs/metric_cache'/audit['cache_key'];target.mkdir(parents=True,exist_ok=True);copied=0;existing=0;metadata=None
 for path in sorted(source.glob('*.npz')):
  with np.load(path,allow_pickle=False) as z:
   meta=str(z['metadata']);r=float(z['r']);h=z['h']
  if metadata is None:metadata=meta
  if meta!=metadata or not np.isfinite(h).all() or h.shape!=(len(theta),4,4):raise ValueError('Invalid auxiliary cache record')
  dest=target/path.name
  if dest.exists():
   with np.load(dest,allow_pickle=False) as z:
    if str(z['metadata'])!=meta or float(z['r'])!=r:raise ValueError('Standard cache has different provenance')
   existing+=1;continue
  temp=dest.with_name(dest.name+'.auxiliary-'+str(os.getpid())+'.tmp');shutil.copyfile(path,temp);os.replace(temp,dest);copied+=1
 missing=[]
 for r in radii:
  path=target/(hashlib.sha256(float(r).hex().encode()).hexdigest()[:24]+'.npz')
  if not path.exists():missing.append(float(r));continue
  with np.load(path,allow_pickle=False) as z:
   if str(z['metadata'])!=metadata or float(z['r'])!=r or z['h'].shape!=(len(theta),4,4) or not np.isfinite(z['h']).all():raise ValueError('Merged cache validation failed')
 record.update(status='all_radius_cache_ready' if not missing else 'tail_cache_merged_waiting_for_primary',audit=audit,copied=copied,already_present=existing,metadata=json.loads(metadata),missing_primary_radii=missing,standard_cache=str(target.relative_to(ROOT)))
 extra.output.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:record[k] for k in ['status','copied','already_present','missing_primary_radii','standard_cache']}),flush=True)
if __name__=='__main__':main()
