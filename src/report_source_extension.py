"""Compare completed source-domain extensions with identical inner samples."""
import argparse
import hashlib
import json
from pathlib import Path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inner',type=Path)
    parser.add_argument('extended',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    old,new=[json.loads(p.read_text()) for p in (args.inner,args.extended)]
    for data in (old,new):
        if data.get('status')!='truncated_single_mode_not_converged' or 'flux' not in data:
            raise ValueError('Both runs must be completed')
        if data.get('flux_validity')=='historical_unreliable_boundary_result':
            raise ValueError('Invalid radial boundary result')
    before=dict(old['parameters']);after=dict(new['parameters'])
    old_panels=before.pop('source_panels');new_panels=after.pop('source_panels')
    if before!=after:raise ValueError('Parameters other than the source domain differ')
    if new_panels[:len(old_panels)]!=old_panels or len(new_panels)<=len(old_panels):
        raise ValueError('New domain must append panels to the old domain')
    count=len(old['samples'])
    if len(new['samples'])<=count:raise ValueError('No added source samples')
    for a,b in zip(old['samples'],new['samples'][:count]):
        if any(a[k]!=b[k] for k in ('r','weight','source')):
            raise ValueError('Inner quadrature samples differ')
        if not b['reused']:raise ValueError('Inner samples were not reused')
    changes={}
    for boundary in ('infinity','horizon'):
        a=old['flux'][boundary]['orbital_energy'];b=new['flux'][boundary]['orbital_energy']
        changes[boundary]=dict(original=a,extended=b,absolute_change=b-a,
                               relative_change=b/a-1 if a else None)
        key='up_coefficient' if boundary=='infinity' else 'z_h'
        z0,z1=complex(*old[key]),complex(*new[key])
        relative=(z1-z0)/z0 if z0 else None
        changes[boundary]['relative_complex_amplitude_change']=(
            [relative.real,relative.imag] if relative is not None else None)
        predicted=(2*relative.real+abs(relative)**2) if relative is not None else None
        changes[boundary]['flux_change_from_amplitude']=predicted
        if a and predicted is not None and abs(predicted-changes[boundary]['relative_change'])>1e-12:
            raise ValueError('Flux and amplitude changes are inconsistent')
    result=dict(status='single_source_extension_not_global_convergence',
        input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (args.inner,args.extended)},
        parameters=before,source_outer_radii=[old_panels[-1],new_panels[-1]],
        reused_inner_samples=count,new_samples=len(new['samples'])-count,flux_changes=changes,
        limitations=['Only the source outer cutoff changes; the Green boundary is fixed',
                     'One extension is not a rigorous bound on all remaining tails',
                     'Added-panel quadrature and all other discretizations still require convergence'])
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(changes,indent=2),flush=True)


if __name__=='__main__':main()
