"""Compare completed batches that differ only in the inner source cutoff."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def compare(original,extended):
    documents=[json.loads(p.read_text()) for p in (original,extended)]
    if any(d['status']!='batch_completed_finite_resolution_not_converged' for d in documents):
        raise ValueError('Both batches must be completed')
    p,q=[d['parameters'] for d in documents]
    if {k:v for k,v in p.items() if k!='source_inner_offset'}!={k:v for k,v in q.items() if k!='source_inner_offset'}:
        raise ValueError('Batch parameters differ beyond source inner cutoff')
    if not 0<q['source_inner_offset']<p['source_inner_offset']:
        raise ValueError('Second cutoff must be smaller and positive')
    maps=[{(c['scalar_ell'],c['scalar_m']):c for c in d['channels']} for d in documents]
    if maps[0].keys()!=maps[1].keys():raise ValueError('Channel sets differ')
    cases=[]
    for mode in sorted(maps[0]):
        paths=[base.parent/channel[mode]['file'] for base,channel in zip((original,extended),maps)]
        data=[json.loads(path.read_text()) for path in paths]
        left,right=[d['parameters'].copy() for d in data]
        panels=[x.pop('source_panels') for x in (left,right)]
        if left!=right or panels[0][1:]!=panels[1][1:] or not panels[1][0]<panels[0][0]:
            raise ValueError('Response parameters or unchanged outer panels differ')
        common=[{s['r']:s for s in d['samples'] if s['r']>panels[0][1]} for d in data]
        if not common[0] or common[0].keys()!=common[1].keys():
            raise ValueError('Unchanged outer quadrature radii differ')
        for radius in common[0]:
            for key in ('weight','source'):
                if common[0][radius][key]!=common[1][radius][key]:
                    raise ValueError('Unchanged outer weight or source differs')
        for d,ch in zip(data,maps):
            if d['flux']!=ch[mode]['flux']:raise ValueError('Batch and response flux differ')
            if not np.isfinite([[s['r'],s['weight'],*s['source']] for s in d['samples']]).all():
                raise ValueError('Nonfinite source samples')
        flux={}
        for boundary in ('horizon','infinity'):
            x,y=[d['flux'][boundary]['orbital_energy'] for d in data]
            flux[boundary]=dict(original=x,extended=y,absolute_change=y-x,relative_change=(y-x)/x if x else None)
        cases.append(dict(ell=mode[0],m=mode[1],inputs=[dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths],
            unchanged_outer_nodes=len(common[0]),flux=flux))
    return dict(status='completed_inner_cutoff_comparison_not_full_convergence',
        inputs=[dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in (original,extended)],
        original_offset=p['source_inner_offset'],extended_offset=q['source_inner_offset'],cases=cases,
        limitation='First logarithmic panel is remeshed. This changes cutoff and its quadrature together; fixed panel order is not an independent quadrature convergence test.')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('original',type=Path);parser.add_argument('extended',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=compare(args.original,args.extended)
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    for c in report['cases']:print(c['ell'],c['m'],c['flux'])


if __name__=='__main__':main()
