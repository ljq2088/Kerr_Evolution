import argparse,hashlib,json
from pathlib import Path
import numpy as np
from lorenz_jet import Jet
from paper_all_component_matching import solve


def main():
    p=argparse.ArgumentParser();p.add_argument('--ellmax',type=int,default=6);p.add_argument('--quadrature',type=int,default=24)
    p.add_argument('--precise-input',action='store_true')
    args=p.parse_args();Jet.coefficient_dtype=np.clongdouble
    if args.precise_input:
        from paper_precise_angular import install_matching_angular
        from paper_analytic_jumps import curvature_jumps
        import paper_all_component_matching
        install_matching_angular()
        paper_all_component_matching.curvature_jumps=curvature_jumps
    d=solve(20.,.8771530275949366,1,args.ellmax,args.quadrature,extended_solve=args.precise_input)
    for k,v in list(d.items()):
        if isinstance(v,np.ndarray):
            if v.dtype==bool:d[k]=v.tolist()
            else:v=np.asarray(v,complex);d[k]=np.stack([v.real,v.imag],axis=-1).tolist()
    root=Path(__file__).resolve().parents[1]
    d['status']='ten_component_local_matching_with_held_out_conditions'
    d['precise_input']=args.precise_input
    d['implementation_sha256']={n:hashlib.sha256((root/'src'/n).read_bytes()).hexdigest() for n in ['paper_all_component_matching.py','paper_full_tetrad.py','paper_precise_angular.py','paper_analytic_jumps.py','paper_jump_basis.py','paper_sourced_matching.py']}
    d['retained_degrees_after_six_mode_guard']=[v for v in d['degree_errors'] if args.ellmax-v['j']>6]
    d['limitations']=['Finite internal ell cutoff; discard the six highest degrees plus the cutoff degree as in the paper',
       'Precise input uses 50-digit angular eigensolve, extended values and analytic curvature jumps; otherwise double inputs. The precise run also solves the same extended-input matrix at 50 digits and records the double QR control',
       'No fitted environmental flux and no author source arrays']
    out=root/'docs/environment_reproduction'/f'all_component_matching_L{args.ellmax}_q{args.quadrature}{"_precise" if args.precise_input else ""}.json'
    out.write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d['degree_errors']),flush=True)


if __name__=='__main__':main()
