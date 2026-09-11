"""Recompute selected actual-source nodes with consistent dense angular factors.

An isolated diagnostic: no cached production tensors or source files are changed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from environment_angular_diagnostic import install_dense_angular_diagnostic
from environment_source import ThresholdCloud,project_source
from environment_lorenz_mode import LorenzMetricMode,ConjugateMetricMode


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs',type=Path,nargs='+')
    parser.add_argument('--sample-indices',type=int,nargs='+',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    reports=[json.loads(path.read_text()) for path in args.inputs]
    common=None
    for data in reports:
        if data.get('status')!='truncated_single_mode_not_converged':
            raise ValueError('Only completed sources can be compared')
        p=data['parameters'];g=p['metric']
        this=(p['alpha'],p['cloud_mass'],g['a'],g['orbital_radius'],abs(g['m_g']),g['ellmax'],p['angular_order'])
        if common is not None and this!=common:raise ValueError('Source parameters differ')
        if p.get('background') is not None or g['m_g']==0:raise ValueError('Requires nonstatic stationary Kerr source')
        common=this
    alpha,mass,a,orbit,mg,ellmax,ntheta=common
    install_dense_angular_diagnostic()
    cloud=ThresholdCloud(alpha=alpha,mass=mass)
    if abs(cloud.a-a)>1e-13:raise ValueError('Cloud spin differs from reference')
    base=LorenzMetricMode(orbit,a,mg,ellmax)
    opposite=ConjugateMetricMode(base)
    rows=[]
    def save():
        result=dict(status='selected_source_nodes_not_flux_convergence',
            parameters=dict(alpha=alpha,cloud_mass=mass,a=a,r0=orbit,abs_metric_m=mg,
                            metric_ellmax=ellmax,angular_order=ntheta,angular_backend='dense-real'),
            inputs=[dict(file=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for path in args.inputs],
            requested_indices=args.sample_indices,completed=rows,
            all_requested_nodes_completed=len(rows)==len(args.inputs)*len(args.sample_indices),
            limitations=['Fresh source amplitudes and homogeneous angular factors use one backend',
                         'Same radial inputs, cloud and quadrature as references',
                         'Selected-node differences do not bound integrated flux errors'])
        temporary=args.output.with_suffix('.tmp')
        temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(args.output)
    save()
    for index in args.sample_indices:
        for path,data in zip(args.inputs,reports):
            p=data['parameters'];sample=data['samples'][index];r=sample['r']
            metric=base if p['metric']['m_g']>0 else opposite
            omega,values=project_source(cloud,[r],orbit,p['scalar_ell'],p['scalar_m'],metric,ntheta=ntheta)
            if abs(omega-p['omega'])>1e-13:raise ValueError('Frequency differs')
            old=complex(*sample['source']);new=values[0]
            row=dict(input=path.name,sample_index=index,r=r,ell=p['scalar_ell'],m=p['scalar_m'],
                original_source=[old.real,old.imag],dense_source=[float(new.real),float(new.imag)],
                absolute_difference=float(abs(new-old)),relative_difference=float(abs(new-old)/abs(old)) if old else None)
            rows.append(row);save();print(row,flush=True)


if __name__=='__main__':main()
