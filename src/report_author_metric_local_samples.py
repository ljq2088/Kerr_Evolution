"""Sample local metric sectors for a direct, unfitted original-author comparison."""
import argparse,json,time
from pathlib import Path
import numpy as np
from pybhpt.swsh import Yslm
from lorenz_metric import spin0_metric,spin1_metric,spin2_metric
from paper_full_tetrad import SPINS
from report_author_hdf5_metric_comparison import weighted_metric,enc,COMPONENTS
from environment_angular_diagnostic import install_dense_angular_diagnostic
from source_provenance import local_dependency_hashes
ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser();p.add_argument('--L',type=int,default=10);p.add_argument('--q',type=int,default=32)
    p.add_argument('--a',type=float,default=.6);p.add_argument('--rp',type=float,default=8.)
    p.add_argument('--m',type=int,default=1);p.add_argument('--output',required=True)
    p.add_argument('--radii',type=float,nargs='+',default=[3.949919837349015,11.986125913940135]);args=p.parse_args()
    out=Path(args.output)
    if out.exists():raise FileExistsError(out)
    install_dense_angular_diagnostic();x,w=np.polynomial.legendre.leggauss(args.q);theta=np.arccos(x)
    angular=np.array([[Yslm(s,j,args.m,theta) if j>=max(abs(s),abs(args.m)) else np.zeros_like(theta) for s in SPINS] for j in range(args.L+1)])
    result=dict(status='sampling',a=args.a,rp=args.rp,m=args.m,lmax=args.L,q=args.q,components=COMPONENTS,
       input_ells=list(range(abs(args.m),args.L+1)),basis='All ten output components spherical; local spin0, spin1, spin2 sectors kept separately',
       implementation_sha256=local_dependency_hashes(ROOT/'src',['report_author_metric_local_samples']),rows=[])
    start=time.perf_counter()
    def save():
        tmp=out.with_suffix('.tmp');tmp.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');tmp.replace(out)
    for r in args.radii:
        samples=np.zeros((3,len(result['input_ells']),args.q,10),complex)
        for k,t in enumerate(theta):
            for j,ell in enumerate(result['input_ells']):
                for sector,fn in [(0,spin0_metric),(1,spin1_metric),(2,spin2_metric)]:
                    if sector==2 and ell<2:continue
                    extra={'full_current':True} if sector==1 else {}
                    g,h=fn(r,t,args.rp,args.a,ell,args.m,order=6,**extra)
                    samples[sector,j,k]=weighted_metric(g,h,r,t,args.a)
            print(f'r={r:.9g},theta={k+1}/{args.q}',flush=True)
        projected=2*np.pi*np.einsum('jck,slkc,k->sljc',angular.conjugate(),samples,w)
        result['rows'].append(dict(r=r,side='In' if r<args.rp else 'Up',local_by_input_ell_and_sector=enc(projected),local_sectors=enc(projected.sum(axis=1)),local=enc(projected.sum(axis=(0,1)))))
        save()
    result.update(status='completed',elapsed_seconds=time.perf_counter()-start);save()
if __name__=='__main__':main()
