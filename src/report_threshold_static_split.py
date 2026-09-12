"""Attribute a saved complete threshold wake to static and nonstatic sectors."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from environment_wake import EnvironmentalWake


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('full_report',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=json.loads(args.full_report.read_text())
    if report['status']!='complete_finite_mode_range_not_paper_reproduction':
        raise ValueError('Require complete finite field data')
    arrays=np.load(args.full_report.with_suffix('.npz'),allow_pickle=False)
    rows=[];fields={}
    root=Path(__file__).resolve().parents[1]
    for dataset in report['datasets']:
        radius=dataset['parameters']['r0'];suffix={41.6:'416',41.8:'418'}[radius]
        selected=[]
        for item in dataset['inputs']:
            path=Path(item['file']);path=path if path.is_absolute() else root/path
            if hashlib.sha256(path.read_bytes()).hexdigest()!=item['sha256']:
                raise ValueError('Input field data changed')
            data=json.loads(path.read_text())
            if data['parameters']['scalar_m']==1:selected.append(data)
        if len(selected)!=5 or sorted(d['parameters']['scalar_ell'] for d in selected)!=[3,5,7,9,11]:
            raise ValueError('Expected five static scalar modes through ell12')
        wake=EnvironmentalWake(selected,ellmax=12,allow_partial=True,allow_mixed_discretization=True)
        edges=arrays['radial_edges_'+suffix];angles=arrays['angular_edges']
        radii=np.sqrt(edges[:-1]*edges[1:])[:,None]
        phi=(angles[:-1]+angles[1:])[None,:]/2
        static=wake.evaluate(radii,np.pi/2,phi)*wake.parameters['alpha']**-3
        full=arrays['field_'+suffix];dynamic=full-static
        if not np.isfinite(static).all():raise ValueError('Nonfinite static field')
        sectors={}
        for name,value in [('full',full),('static_m1',static),('nonstatic',dynamic)]:
            index=np.unravel_index(np.argmax(abs(value)),value.shape)
            sectors[name]=dict(maximum=float(np.max(abs(value))),
                maximum_r=float(radii[index[0],0]),maximum_phi=float(phi[0,index[1]]),
                rms=float(np.sqrt(np.mean(abs(value)**2))))
            fields[name+'_'+suffix]=value
        rows.append(dict(orbit=radius,sectors=sectors))
        print(rows[-1],flush=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(args.output.with_suffix('.npz'),**fields)
    args.output.write_text(json.dumps(dict(status='sector_attribution_not_error_estimate',
        full_report=str(args.full_report),sha256=hashlib.sha256(args.full_report.read_bytes()).hexdigest(),
        rows=rows,limitation='Nonstatic field is the complex full field minus the static field, before taking absolute values. Static terms are retained in the physical result. RMS uses plotting cells, not a physical volume measure.'),indent=2)+'\n')


if __name__=='__main__':main()
