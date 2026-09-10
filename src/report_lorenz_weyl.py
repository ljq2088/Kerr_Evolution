"""Audit raw and time-integrated Weyl conventions against the public table."""
import json
import tarfile
import hashlib
from pathlib import Path
from importlib.metadata import version
from lorenz_weyl import weyl_amplitudes


def main():
    root=Path(__file__).resolve().parents[1]
    archive=root/'outputs/lorenz_reference/2406.12510v3.tar'
    with tarfile.open(archive) as tar:
        lines=tar.extractfile('amplitudes.dat').read().decode().splitlines()
    header=lines[0].split()
    cases=[]
    for line in lines[1:]:
        values=line.split()
        r=float(values[0])
        if r not in (4.,6.,10.,20.):
            continue
        data=dict(zip(header,values))
        raw=weyl_amplitudes(r)
        omega=2/(r**1.5+.6)
        for s,names in ((-2,('Psi4Inf','Psi4Hor')),(2,('Psi0Inf','Psi0Hor'))):
            ref=[complex(data[k].replace('*I','j')) for k in names]
            integrated=[1j*z/omega for z in raw[s]]
            cases.append(dict(r0=r,s=s,omega=omega,
                              raw_curvature=[[z.real,z.imag] for z in raw[s]],
                              integrated_curvature=[[z.real,z.imag] for z in integrated],
                              reference=[[z.real,z.imag] for z in ref],
                              raw_relative_error=[abs(x-y)/abs(y) for x,y in zip(raw[s],ref)],
                              integrated_relative_error=[abs(x-y)/abs(y) for x,y in zip(integrated,ref)]))
    report=dict(status='Weyl_inputs_only_not_full_Lorenz_metric',pybhpt_version=version('pybhpt'),
                source='https://arxiv.org/src/2406.12510v3',member='amplitudes.dat',
                archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
                convention='Table matches inverse time derivative i/omega of raw pybhpt Weyl amplitudes; retain both.',
                cases=cases)
    (root/'docs/environment_reproduction/lorenz_weyl_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Maximum relative integrated-amplitude error:',
          max(max(c['integrated_relative_error']) for c in cases))


if __name__=='__main__':
    main()
