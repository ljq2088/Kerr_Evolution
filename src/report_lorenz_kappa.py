"""Independent step, tolerance and outer-boundary audit of resolvent kappa."""
import json
from pathlib import Path
import numpy as np
from lorenz_ghp import KerrGHP
from lorenz_kappa import kappa_jet


def main():
    configurations=[(1e-4,1e-12,2000.),(5e-5,1e-12,2000.),
                    (2.5e-5,1e-12,2000.),(5e-5,2e-13,2000.),
                    (5e-5,1e-12,1000.)]
    cases=[]
    for r in (4.5,8.):
        g=KerrGHP(r,1.1,.6,omega=2/(6**1.5+.6),m=2,order=4)
        values=[]
        for step,rtol,rmax in configurations:
            k=kappa_jet(g,6.,2,step,rtol,rmax)
            values.append(np.array([k.value,k.derivative(0).value,k.derivative(1).value]))
        reference=values[2]
        for config,value in zip(configurations,values):
            case=dict(r=r,step=config[0],rtol=config[1],rmax=config[2],
                      kappa=[float(value[0].real),float(value[0].imag)],
                      relative_value_and_gradient_change=float(np.max(abs(value-reference))/np.max(abs(reference))))
            cases.append(case)
            print(case,flush=True)
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/kappa_convergence.json'
    out.write_text(json.dumps(dict(status='vacuum_resolvent_only',reference='step=2.5e-5,rtol=1e-12,rmax=2000',cases=cases),indent=2)+'\n')


if __name__=='__main__':
    main()
