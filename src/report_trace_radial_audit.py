"""Compare independent scalar resolvent and Teukolsky radial trace solutions."""
import json
from pathlib import Path
import numpy as np
from pybhpt.radial import RadialTeukolsky
from lorenz_kappa import scalar_resolvent


def main():
    cases=[]
    for ell in (2,4,6,8,10):
        radial,zi,zh=scalar_resolvent(6.,.6,ell,2,0.)
        radii=np.array([4.5,6.,8.])
        mode=RadialTeukolsky(0,ell,2,.6,2/(6**1.5+.6),radii)
        mode.solve()
        ri=np.array([mode.radialsolution('In',i) for i in range(3)])
        ru=np.array([mode.radialsolution('Up',i) for i in range(3)])
        di=mode.radialderivative('In',1)
        du=mode.radialderivative('Up',1)
        w=(36-12+.36)*(ri[1]*du-ru[1]*di)
        independent=np.array([ri[0]*ru[1],ri[1]*ru[1],ri[1]*ru[2]])/w
        direct=np.array([radial.insol.sol(4.5)[0]*radial.upsol.sol(6.)[0],
                         radial.insol.sol(6.)[0]*radial.upsol.sol(6.)[0],
                         radial.insol.sol(6.)[0]*radial.upsol.sol(8.)[0]])/radial.w0
        case=dict(ell=ell,relative_green_difference=float(np.max(abs(direct/independent-1))),
                  relative_wronskian_variation=float(np.max(abs(radial.wronskian(radii)/radial.w0-1))))
        cases.append(case)
        print(case,flush=True)
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/trace_radial_audit.json'
    out.write_text(json.dumps(dict(cases=cases),indent=2)+'\n')


if __name__=='__main__':
    main()
