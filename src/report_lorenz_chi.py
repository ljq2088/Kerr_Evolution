"""Compact auxiliary scalar versus the public table; kappa is separate."""
import json
import tarfile
from pathlib import Path
from lorenz_chi import chi_amplitudes


def main():
    root=Path(__file__).resolve().parents[1]
    with tarfile.open(root/'outputs/lorenz_reference/2406.12510v3.tar') as tar:
        lines=tar.extractfile('amplitudes.dat').read().decode().splitlines()
    keys=lines[0].split()
    cases=[]
    for line in lines[1:]:
        vals=line.split()
        r=float(vals[0])
        if r not in (4.,6.,10.,20.):
            continue
        row=dict(zip(keys,vals))
        ref=[complex(row[k].replace('*I','j')) for k in ('ChiInf','ChiHor')]
        omega=2/(r**1.5+.6)
        raw=chi_amplitudes(r)
        mapped=[-1j*z/(2*omega) for z in raw]
        cases.append(dict(r0=r,raw=[[z.real,z.imag] for z in raw],
                          table_convention=[[z.real,z.imag] for z in mapped],
                          reference=[[z.real,z.imag] for z in ref],
                          mapped_relative_error=[abs(x-y)/abs(y) for x,y in zip(mapped,ref)]))
    report=dict(status='compact_source_chi_only_kappa_not_included',
                source='https://arxiv.org/src/2406.12510v3',member='amplitudes.dat',
                table_mapping='-i/(2 omega) times raw compact scalar',cases=cases)
    (root/'docs/environment_reproduction/lorenz_chi_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Maximum table-mapped relative error:',max(max(c['mapped_relative_error']) for c in cases))


if __name__=='__main__':
    main()
