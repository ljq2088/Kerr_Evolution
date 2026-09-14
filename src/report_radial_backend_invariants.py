"""Compare normalization-invariant Teukolsky Green kernels across backends."""
import json
from pathlib import Path
import numpy as np
from pybhpt.radial import RadialTeukolsky,available_methods


def main():
    a=.8771530275949366;r0=20.;m=1;omega=m/(r0**1.5+a)
    radii=np.array([2.,10.,20.,20.00005,40.]);rows=[]
    folder=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    print('Available',available_methods(),flush=True)
    for s in (-2,-1,1,2):
      for ell in (max(abs(s),1),4,8,14):
        baseline=None
        for method,tolerance in (('AUTO',None),('AUTO',1e-13),('MST' if abs(s)==2 else 'HBL',1e-13),('TEUK',1e-13)):
          radial=RadialTeukolsky(s,ell,m,a,omega,radii);radial.solve(method=method,rtol=tolerance)
          ri=np.array([radial.radialsolution('In',i) for i in range(len(radii))]);di=np.array([radial.radialderivative('In',i) for i in range(len(radii))])
          ru=np.array([radial.radialsolution('Up',i) for i in range(len(radii))]);du=np.array([radial.radialderivative('Up',i) for i in range(len(radii))])
          W=(r0*r0-2*r0+a*a)**(s+1)*(ri[2]*du[2]-ru[2]*di[2])
          kernel=np.where(radii<r0,ru[2]*ri,ri[2]*ru)/W
          logder=np.concatenate([di/ri,du/ru])
          if not np.all(np.isfinite(kernel)) or not np.all(np.isfinite(logder)):raise ValueError('Nonfinite backend result')
          if baseline is None:baseline=(kernel,logder)
          row=dict(spin=s,ell=ell,method=method,rtol=tolerance,green_kernel_max_relative_change=float(np.max(abs(kernel/baseline[0]-1))),
                   logarithmic_derivative_max_relative_change=float(np.max(abs(logder/baseline[1]-1))))
          rows.append(row);print(json.dumps(row),flush=True)
          result=dict(status='radial_backend_invariants_diagnostic_not_flux_validation',a=a,r0=r0,m=m,omega=omega,radii=radii.tolist(),rows=rows,
            limitations=['No amplitude normalization comparison; uses invariant combinations',
                         'Source differential operators and reconstructed metric not tested here',
                         'MST restricted to spin2 after preliminary spin-1 ell1 returned nonfinite; spin1 uses HBL'])
          out=folder/'radial_backend_invariants.json';tmp=out.with_suffix('.tmp');tmp.write_text(json.dumps(result,indent=2)+'\n');tmp.replace(out)


if __name__=='__main__':main()
