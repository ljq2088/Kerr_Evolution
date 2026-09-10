"""Resumable actual-source quadrature for one nonstatic scalar response mode.

This deliberately records source cutoffs and finite angular truncation. No
output here is a converged paper flux until those limits have been varied.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from environment_source import ThresholdCloud,project_source
from environment_lorenz_mode import LorenzMetricMode
from environment_radial import RadialGreen
from environment_cloud import mode_flux


def encode(z):
    return [float(z.real),float(z.imag)]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--radial-order',type=int,default=2)
    parser.add_argument('--angular-order',type=int,default=6)
    parser.add_argument('--metric-ellmax',type=int,default=4)
    parser.add_argument('--source-inner-offset',type=float,default=.05)
    parser.add_argument('--source-outer-radius',type=float,default=320.)
    args=parser.parse_args()
    if min(args.radial_order,args.angular_order)<2:
        raise ValueError('Quadrature orders must be at least two')
    cloud=ThresholdCloud(alpha=.3)
    metric=LorenzMetricMode(20.,cloud.a,2,args.metric_ellmax)
    omega=cloud.omega+metric.omega
    green=RadialGreen(cloud.a,cloud.mu,omega,3,3,rmax=1000.,offset=1e-4,rtol=1e-11)
    inner=cloud.rp+args.source_inner_offset
    outer=args.source_outer_radius
    if not green.rmin<=inner<20.<outer<=min(green.rmax,cloud.rmax):
        raise ValueError('Source cutoffs must enclose the orbit and lie inside both solved domains')
    panels=np.array([inner]+[r for r in (3.,6.,12.,20.,40.,80.,160.,320.) if inner<r<outer]+[outer])
    x,w=np.polynomial.legendre.leggauss(args.radial_order)
    radii=np.concatenate([(a+b)/2+(b-a)*x/2 for a,b in zip(panels[:-1],panels[1:])])
    weights=np.concatenate([(b-a)*w/2 for a,b in zip(panels[:-1],panels[1:])])
    directory=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    suffix=''
    if args.source_inner_offset!=.05 or outer!=320.:
        suffix=f'_inner{args.source_inner_offset:g}_outer{outer:g}'
    out=directory/f'forced_mode_nr{args.radial_order}_nt{args.angular_order}_L{args.metric_ellmax}{suffix}.json'
    metadata=dict(alpha=.3,cloud_mass=1.,metric=metric.provenance,scalar_ell=3,scalar_m=3,
                  omega=float(omega),radial_order=args.radial_order,angular_order=args.angular_order,
                  source_panels=panels.tolist(),green_outer_radius=1000.,green_horizon_offset=1e-4)
    samples=[]
    if out.exists():
        previous=json.loads(out.read_text())
        if previous['parameters']!=metadata:
            raise ValueError('Saved source samples have different physical or numerical parameters')
        samples=previous['samples']
        if len(samples)>len(radii) or any(row['r']!=float(radii[i]) for i,row in enumerate(samples)):
            raise ValueError('Saved quadrature nodes differ')
    def save():
        temporary=out.with_suffix('.tmp')
        temporary.write_text(json.dumps(result,indent=2)+'\n')
        temporary.replace(out)
    result=dict(status='source_sampling_in_progress',parameters=metadata,samples=samples)
    for index in range(len(samples),len(radii)):
        r=float(radii[index])
        frequency,J=project_source(cloud,[r],20.,3,3,metric,ntheta=args.angular_order)
        if abs(frequency-omega)>1e-14:
            raise RuntimeError('Source and Green frequencies differ')
        samples.append(dict(r=r,weight=float(weights[index]),source=encode(J[0])))
        save()
        print(f'Source radius {index+1}/{len(radii)}: r={r:.8g}',flush=True)
    J=np.array([complex(*sample['source']) for sample in samples])
    ingoing=green.insol.sol(radii)[0]
    outgoing=green.upsol.sol(radii)[0]
    zi=np.sum(weights*ingoing*J)/green.w0
    zh=np.sum(weights*outgoing*J)/green.w0
    field_samples=[]
    cutoff_sequence=[]
    for cutoff in panels[1:]:
        selected=radii<cutoff
        ci=np.sum(weights[selected]*ingoing[selected]*J[selected])/green.w0
        ch=np.sum(weights[selected]*outgoing[selected]*J[selected])/green.w0
        cutoff_sequence.append(dict(outer_source_cutoff=float(cutoff),z_inf=encode(ci),z_h=encode(ch)))
    for r in np.concatenate(([green.rmin],panels,[green.rmax])):
        left=radii<r
        cu=np.sum(weights[left]*ingoing[left]*J[left])/green.w0
        cv=np.sum(weights[~left]*outgoing[~left]*J[~left])/green.w0
        u,du=green.insol.sol(r)
        v,dv=green.upsol.sol(r)
        field_samples.append(dict(r=float(r),field=encode(v*cu+u*cv),derivative=encode(dv*cu+du*cv)))
    result.update(status='truncated_single_mode_not_converged',z_inf=encode(zi),z_h=encode(zh),
                  radial_response_at_panel_boundaries=field_samples,
                  outer_source_cutoff_sequence=cutoff_sequence,
                  flux_scaling='per q^2*(cloud mass/M); unit Killing cloud mass, not alpha^-3 rescaled',
                  flux=mode_flux(omega,3,cloud.omega,cloud.m,cloud.mu,cloud.a,zi,zh),
                  wronskian_relative_spread=float(np.max(abs(green.wronskian(radii)/green.w0-1))),
                  missing_convergence=['radial quadrature','horizon source cutoff','outer source cutoff',
                                       'angular projection','metric ell truncation','remaining m modes and static completion'])
    save()
    print(json.dumps({k:result[k] for k in ('status','z_inf','z_h','flux','wronskian_relative_spread')}),flush=True)


if __name__=='__main__':
    main()
