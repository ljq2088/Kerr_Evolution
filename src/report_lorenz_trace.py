"""Compare independent trace calculation with public 2406.12510v3 data."""
import json
import hashlib
import tarfile
from pathlib import Path
from lorenz_trace import trace_amplitudes


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
        reference=[complex(data[k].replace('*I','j')) for k in ('hInf','hHor')]
        baseline=trace_amplitudes(r,rmax=2000.,offset=1e-5)
        refined=trace_amplitudes(r,rmax=4000.,offset=1e-6)
        cases.append(dict(r0=r,reference=[[z.real,z.imag] for z in reference],
                          computed=[[z.real,z.imag] for z in refined],
                          absolute_errors=[abs(x-y) for x,y in zip(refined,reference)],
                          relative_errors=[abs(x-y)/abs(y) for x,y in zip(refined,reference)],
                          cutoff_change=[abs(x-y) for x,y in zip(baseline,refined)]))
    report=dict(source='https://arxiv.org/src/2406.12510v3',member='amplitudes.dat',
                archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
                a=.6,ell=2,m=2,status='trace_only_not_full_Lorenz_metric',cases=cases)
    out=root/'docs/environment_reproduction/lorenz_trace_validation.json'
    out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
