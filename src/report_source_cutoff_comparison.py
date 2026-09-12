"""Compare identical inner source samples with an appended outer source panel."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from environment_response import SampledResponse


def compare(old,new,radii):
    for data in (old,new):
        if data.get('status')!='truncated_single_mode_not_converged':
            raise ValueError('Both source integrations must be complete')
    a=dict(old['parameters']);b=dict(new['parameters'])
    pa=a.pop('source_panels');pb=b.pop('source_panels')
    if a!=b:raise ValueError('Other physical or numerical parameters differ')
    if pb[:len(pa)]!=pa or len(pb)<=len(pa):
        raise ValueError('New source panels must strictly append the old panels')
    lookup={row['r']:row for row in new['samples']}
    if len(lookup)!=len(new['samples']):raise ValueError('Duplicate source radius')
    for row in old['samples']:
        other=lookup.get(row['r'])
        if other is None:raise ValueError('Missing original source radius')
        for key in ('r','weight','source'):
            if row[key]!=other[key]:raise ValueError(f'Original sample {key} differs')
    if np.max(radii)>pa[-1]:raise ValueError('Comparison must stay inside the original source cutoff')
    fields=[SampledResponse.from_report(data).evaluate(radii)[0] for data in (old,new)]
    if not all(np.all(np.isfinite(v)) for v in fields):raise ValueError('Nonfinite comparison field')
    scale=float(np.max(abs(fields[0])))
    if scale==0:raise ValueError('Cannot normalize a vanishing field')
    fluxes={}
    for boundary in ('infinity','horizon'):
        x,y=[data['flux'][boundary]['orbital_energy'] for data in (old,new)]
        fluxes[boundary]=dict(original=x,extended=y,absolute_change=y-x,
                              relative_change=y/x-1 if x else None)
    return dict(ell=a['scalar_ell'],m=a['scalar_m'],
        original_nodes=len(old['samples']),extended_nodes=len(new['samples']),
        original_cutoff=pa[-1],extended_cutoff=pb[-1],
        exact_inner_source_reuse=True,fluxes=fluxes,
        field_relative_sup_change=float(np.max(abs(fields[1]-fields[0]))/scale))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('original_batch',type=Path)
    parser.add_argument('extended_batch',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--radial-points',type=int,default=160)
    args=parser.parse_args()
    batches=[json.loads(path.read_text()) for path in (args.original_batch,args.extended_batch)]
    if any(b['status']!='batch_completed_finite_resolution_not_converged' for b in batches):
        raise ValueError('Both batches must be complete')
    maps=[{(c['scalar_ell'],c['scalar_m']):c['file'] for c in b['channels']} for b in batches]
    if maps[0].keys()!=maps[1].keys():raise ValueError('Batch modes differ')
    rows=[];radii=None
    for key in sorted(maps[0]):
        paths=[batch.parent/maps[i][key]
               for i,batch in enumerate((args.original_batch,args.extended_batch))]
        data=[json.loads(path.read_text()) for path in paths]
        if radii is None:
            p=data[0]['parameters'];rp=1+np.sqrt(1-p['metric']['a']**2)
            if args.radial_points<2:raise ValueError('At least two comparison radii required')
            radii=np.geomspace(max(rp+.06,p['source_panels'][0]),p['source_panels'][-1],args.radial_points)
        row=compare(*data,radii)
        row['inputs']={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
        rows.append(row)
        print(key,row['field_relative_sup_change'],flush=True)
    report=dict(status='source_cutoff_comparison_not_global_convergence',
        field_radii=radii.tolist(),rows=rows,
        field_metric='max abs difference / max abs original field at listed radii',
        limitation='No interpolation-order, metric-truncation, all-mode or paper-agreement claim.')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
