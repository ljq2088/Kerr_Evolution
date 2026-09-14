"""Compare direct sourced 2023 kappa with the mapped 2024 resolvent."""
import hashlib,json
from pathlib import Path
import numpy as np
from paper_kappa_radial import PaperKappaRadial
from report_paper_kappa_tables import calculate
from environment_trace_variation import TraceMassVariation
from lorenz_chi import chi_amplitudes
from lorenz_metric import _homogeneous_radial_data
from environment_angular_diagnostic import install_dense_angular_diagnostic


def enc(value):
    value=np.asarray(value,complex)
    return np.stack([value.real,value.imag],axis=-1).tolist()


def compare(a,r0,m,ell,rmax,radii):
    w=m/(r0**1.5+a);j0,j1,_=calculate(a,r0,m,ell)
    direct=PaperKappaRadial(r0,a,m,ell,[j0,j1],rmax=rmax)
    variation=TraceMassVariation(r0,a,ell,m,rmax=8000.)
    rows=[]
    for r in radii:
        f,fp,df,dfp=variation.radial_state(r)
        bc,index=('In',1) if r<r0 else ('Up',0)
        _,R,Rp=_homogeneous_radial_data(0,ell,m,a,w,r,bc)
        chi=chi_amplitudes(r0,a,ell,m)[index]*np.array([R,Rp])
        ref=np.array([df,dfp])/2-1j*chi/(2*w)
        result=direct.state(r)
        rows.append(dict(r=r,direct_trace=enc(result[:2]),mapped_trace=enc([f,fp]),
            direct_kappa=enc(result[2:]),mapped_kappa=enc(ref),
            trace_relative_max=float(np.max(abs(result[:2]/[f,fp]-1))),
            kappa_relative_max=float(np.max(abs(result[2:]/ref-1)))))
    return dict(a=a,r0=r0,m=m,ell=ell,rmax=rmax,jumps=enc([j0,j1]),
                maximum_relative_kappa_difference=max(x['kappa_relative_max'] for x in rows),rows=rows)


def main():
    install_dense_angular_diagnostic();rows=[]
    cases=[(0.,6.,1,1,4000.,[2.001,3.,5.9,6.1,10.,40.]),
           (.6,6.,2,2,4000.,[1.801,3.,5.9,6.1,10.,40.])]
    cases += [(.8771530275949366,20.,1,1,b,[1.481,3.,10.,19.9,20.1,40.,200.]) for b in (2000.,4000.,8000.)]
    for case in cases:
        row=compare(*case);rows.append(row)
        print(row['a'],row['rmax'],row['maximum_relative_kappa_difference'],flush=True)
    root=Path(__file__).resolve().parents[1]
    sources=['src/paper_kappa_radial.py','src/report_paper_kappa_radial.py','src/environment_trace_variation.py']
    result=dict(status='direct_2023_kappa_boundary_problem_compared_to_2024_resolvent',rows=rows,
      input_sha256={f:hashlib.sha256((root/f).read_bytes()).hexdigest() for f in sources},
      limitations=['Same physical trace source and prescribed jumps; independent endpoint series and coupled ODE',
                   'Tests diagonal particular solution plus matched homogeneous compensation',
                   'No total environmental flux change or author radial data comparison'])
    (root/'docs/environment_reproduction/paper_kappa_radial_audit.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
