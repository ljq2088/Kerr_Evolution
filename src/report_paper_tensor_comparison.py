"""Compare three explicit 2023 tensor components with the production spin-2 sector."""
import json
from pathlib import Path
import numpy as np
from lorenz_metric import spin2_metric,_homogeneous_radial_data
from lorenz_weyl import weyl_amplitudes
from paper_sourced_matching import tensor_components
from report_nonstatic_paper_projection import transform,PAIRS
from environment_angular_diagnostic import install_dense_angular_diagnostic


def main():
    install_dense_angular_diagnostic();rows=[]
    for a,r0,m,ell,r in [(0.,6.,2,2,4.),(.6,6.,2,2,4.)]+[(.8771530275949366,20.,1,2,r) for r in (3.,10.,40.)]:
        theta=1.1;g,h=spin2_metric(r,theta,r0,a,ell,m,10)
        bc,index=('In',1) if r<r0 else ('Up',0)
        _,R,Rp=_homogeneous_radial_data(-2,ell,m,a,g.omega,r,bc)
        amp=-4*weyl_amplitudes(r0,a,ell,m)[-2][index]
        direct=np.array([v.value for v in tensor_components(g,ell,amp*R,amp*Rp)])
        T,_=transform(r,theta,a)
        reference=(T@np.array([h[i][j].value for i,j in PAIRS]))[[0,2,4]]
        enc=lambda z:np.stack([z.real,z.imag],axis=-1).tolist()
        rows.append(dict(a=a,r0=r0,m=m,ell=ell,r=r,theta=theta,explicit=enc(direct),
            production=enc(reference),relative_max=float(np.max(abs(direct-reference))/np.max(abs(reference)))))
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/paper_tensor_comparison.json'
    out.write_text(json.dumps(dict(status='three_tensor_components_not_full_metric_validation',rows=rows,
        limitations=['Shared curvature radial data, independent reconstruction formulas',
                     'Double precision, only the specified modes and points']),indent=2)+'\n')
    print(json.dumps(rows,indent=2))


if __name__=='__main__':main()
