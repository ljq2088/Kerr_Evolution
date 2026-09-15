"""Compare local metric sectors against unmodified public author HDF5 data.

Reference is the ConorDyson/KerrLorenzMSF 2026 fork sample at a=.6,rp=8,
not the unpublished arrays used for the 2025 environmental-flux figures.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
import h5py
from pybhpt.swsh import Yslm
from lorenz_metric import spin0_metric,spin1_metric,spin2_metric
from paper_full_tetrad import project_metric,SPINS
from environment_angular_diagnostic import install_dense_angular_diagnostic
from source_provenance import local_dependency_hashes
ROOT=Path(__file__).resolve().parents[1]
COMPONENTS=['h_l+l+','h_l-l-','h_m+m+','h_m-m-','rho_h_l+m+','rhob_h_l+m-','rhob_h_l-m+','rho_h_l-m-','sigma_delta_h_l+l-','h']

def enc(z):
    v=np.asarray(z,complex);return np.stack([v.real,v.imag],axis=-1).tolist()

def weighted_metric(g,h,r,t,a):
    H=np.array([[x.value for x in row] for row in h],complex)
    v=project_metric(H,r,t,a);rho=r+1j*a*np.cos(t);rc=rho.conjugate()
    sigma=r*r+a*a*np.cos(t)**2;delta=r*r-2*r+a*a
    return np.array([*v[:4],rho*v[4],rc*v[5],rc*v[6],rho*v[7],sigma*delta*v[8],(delta*v[8]+v[9])/sigma])

def errors(value,reference):
    out=[]
    for c,name in enumerate(COMPONENTS):
        v=value[:,c];ref=reference[:,c];scale=np.linalg.norm(ref)
        out.append(dict(component=name,reference_l2=float(scale),difference_l2=float(np.linalg.norm(v-ref)),
                        relative_l2=float(np.linalg.norm(v-ref)/scale) if scale>1e-20 else None,
                        largest_reference_coefficient_ratio=enc(v[np.argmax(abs(ref))]/ref[np.argmax(abs(ref))]) if scale>1e-20 else None))
    return out

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--m',type=int,default=1)
    parser.add_argument('--q',type=int,default=24);parser.add_argument('--radii',type=float,nargs='+',default=[4.,12.])
    parser.add_argument('--output',required=True);args=parser.parse_args()
    out=Path(args.output)
    if out.exists():raise FileExistsError('Use fresh output')
    source=ROOT/f'outputs/paper_metric_reference/ConorDyson_KerrLorenzMSF/GenerationCodes/Numerics-Asymptotics/h1-mmode-gen/h1mmodedat/data600/data600-80/data/h1_a0.6_rp8.0_l4_m{args.m}.h5'
    install_dense_angular_diagnostic();x,w=np.polynomial.legendre.leggauss(args.q);theta=np.arccos(x)
    angular=np.array([[Yslm(s,j,args.m,theta) if j>=max(abs(s),abs(args.m)) else np.zeros_like(theta) for s in SPINS] for j in range(5)])
    result=dict(status='sampling',a=.6,rp=8.,m=args.m,lmax=4,q=args.q,components=COMPONENTS,
                reference_path=str(source),reference_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                reference_repository='https://github.com/ConorDyson/KerrLorenzMSF',
                reference_manifest_sha256=hashlib.sha256((ROOT/'outputs/paper_metric_reference/manifest.json').read_bytes()).hexdigest(),
                implementation_sha256=local_dependency_hashes(ROOT/'src',['report_author_hdf5_metric_comparison']),
                scope='Author sample metric data at a=.6,rp8,L4; neither actual paper flux parameters nor a converged-high-L reference',
                limitations=['HDF5 contains no code commit/physical parameter attributes; parameters checked against repository filenames/configurations',
                             'Current fork contains inconsistent generation-code copies; HDF5 generation commit is unknown',
                             'No fitted scale, phase, or mode substitution is applied'],rows=[])
    start=time.perf_counter()
    def save():
        tmp=out.with_suffix('.tmp');tmp.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');tmp.replace(out)
    with h5py.File(source) as f:
        group=f[f'm_{args.m}']
        for target in args.radii:
            side='In' if target<8 else 'Up';grid=group['r_in' if side=='In' else 'r_up'][()]
            index=int(np.argmin(abs(grid-target)));r=float(grid[index])
            if r==8:raise ValueError('Use an off-orbit sample for this comparison')
            ref=np.array([group[side][c][:,index]['Re']+1j*group[side][c][:,index]['Im'] for c in COMPONENTS]).T
            samples=np.zeros((3,args.q,10),complex);dipole=np.zeros((args.q,10),complex)
            for k,t in enumerate(theta):
                for ell in range(abs(args.m),5):
                    for sector,fn in [(0,spin0_metric),(1,spin1_metric),(2,spin2_metric)]:
                        if sector==2 and ell<2:continue
                        extra={'full_current':True} if sector==1 else {}
                        g,h=fn(r,t,8.,.6,ell,args.m,order=6,**extra)
                        value=weighted_metric(g,h,r,t,.6);samples[sector,k]+=value
                        if ell==1 and sector==1:dipole[k]+=value
                print(f'm={args.m},r={r:.9g},theta={k+1}/{args.q}',flush=True)
            projected=2*np.pi*np.einsum('jck,skc,k->sjc',angular.conjugate(),samples,w)
            dip=2*np.pi*np.einsum('jck,kc,k->jc',angular.conjugate(),dipole,w)
            total=np.sum(projected,axis=0)
            row=dict(r=r,side=side,reference_grid_index=index,reference=enc(ref),local=enc(total),local_sectors=enc(projected),
                     components_error=errors(total,ref))
            if args.m==1:row.update(local_omitting_dipole_vector=enc(total-dip),components_error_omitting_dipole_vector=errors(total-dip,ref))
            result['rows'].append(row);save()
            print('component relative errors',[v['relative_l2'] for v in row['components_error']],flush=True)
    result.update(status='completed_author_sample_comparison',elapsed_seconds=time.perf_counter()-start);save()

if __name__=='__main__':main()
