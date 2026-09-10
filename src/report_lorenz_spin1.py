"""Convention audit; keeps raw results and table-normalized values explicit."""
import json
import tarfile
from pathlib import Path
from lorenz_spin1 import spin1_amplitudes


def main():
    root=Path(__file__).resolve().parents[1]
    with tarfile.open(root/'outputs/lorenz_reference/2406.12510v3.tar') as tar:
        lines=tar.extractfile('amplitudes.dat').read().decode().splitlines()
    header=lines[0].split()
    cases=[]
    for line in lines[1:]:
        cols=line.split()
        r=float(cols[0])
        if r not in (4.,6.,10.,20.):
            continue
        data=dict(zip(header,cols))
        omega=2/(r**1.5+.6)
        for spin,names in ((-1,('Phi2Inf','Phi2Hor')),(1,('Phi0Inf','Phi0Hor'))):
            ref=[complex(data[k].replace('*I','j')) for k in names]
            raw=spin1_amplitudes(r,spin=spin)
            ratios=[expected/computed for expected,computed in zip(ref,raw)]
            mapped=[2**.5*z/omega**2 for z in raw]
            cases.append(dict(r0=r,spin=spin,omega=omega,
                              computed=[[z.real,z.imag] for z in raw],reference=[[z.real,z.imag] for z in ref],
                              reference_over_computed=[[z.real,z.imag] for z in ratios],
                              table_convention=[[z.real,z.imag] for z in mapped],
                              mapped_relative_error=[abs(x-y)/abs(y) for x,y in zip(mapped,ref)],
                              ratio_times_omega_squared=[[z.real*omega**2,z.imag*omega**2] for z in ratios]))
    out=root/'docs/environment_reproduction/lorenz_spin1_validation.json'
    out.write_text(json.dumps(dict(status='source_amplitudes_and_table_bridge_only',
                                  source='https://arxiv.org/src/2406.12510v3',member='amplitudes.dat',
                                  table_mapping='sqrt(2)/omega^2 times raw spin-1 amplitudes',
                                  cases=cases),indent=2)+'\n')
    for c in cases:
        print(c['r0'],c['spin'],max(c['mapped_relative_error']))


if __name__=='__main__':
    main()
