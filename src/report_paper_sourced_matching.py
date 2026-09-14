"""Reproducible CLI for the diagnostic 2023 particle-source matching solve."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from lorenz_jet import Jet
from paper_sourced_matching import solve_jumps


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--a',type=float,required=True)
    parser.add_argument('--r0',type=float,required=True)
    parser.add_argument('--m',type=int,required=True)
    parser.add_argument('--ellmax',type=int,required=True)
    parser.add_argument('--quadrature',type=int,default=20)
    parser.add_argument('--order',type=int,default=10)
    parser.add_argument('--extended',action='store_true')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.extended:Jet.coefficient_dtype=np.clongdouble
    result=solve_jumps(args.r0,args.a,args.m,args.ellmax,args.quadrature,args.order)
    for key,value in list(result.items()):
        if isinstance(value,np.ndarray):
            value=np.asarray(value,complex)
            result[key]=np.stack([value.real,value.imag],axis=-1).tolist()
    root=Path(__file__).resolve().parents[1]
    result['jet_precision']='extended coefficients, double eigenfunctions/radial inputs/QR' if args.extended else 'double'
    result['status']='local_three_component_matching_diagnostic_not_full_metric'
    result['input_sha256']={f:hashlib.sha256((root/f).read_bytes()).hexdigest() for f in
        ['src/paper_jump_basis.py','src/paper_sourced_matching.py']}
    result['limitations']=['Highest retained angular modes can fail unused matching equations',
        'Not an arbitrary-precision solve; do not infer accuracy from selected-equation residual',
        'No flux normalization or fitted paper flux enters this solve']
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['condition_number','selected_residual','unused_residual']}))


if __name__=='__main__':main()
