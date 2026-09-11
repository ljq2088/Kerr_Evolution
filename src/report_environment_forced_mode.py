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
    parser.add_argument('--alpha',type=float,default=.3)
    parser.add_argument('--background',choices=('threshold-kerr','schwarzschild-frozen'),default='threshold-kerr')
    parser.add_argument('--orbital-radius',type=float,default=20.)
    parser.add_argument('--green-horizon-offset',type=float,default=1e-4)
    parser.add_argument('--green-outer-radius',type=float,default=1000.)
    parser.add_argument('--infinity-method',choices=('series','coulomb'),default='series')
    parser.add_argument('--radial-order',type=int,default=2)
    parser.add_argument('--angular-order',type=int,default=6)
    parser.add_argument('--metric-ellmax',type=int,default=4)
    parser.add_argument('--metric-m',type=int,default=2)
    parser.add_argument('--scalar-ell',type=int)
    parser.add_argument('--source-inner-offset',type=float,default=.05)
    parser.add_argument('--source-outer-radius',type=float,default=320.)
    parser.add_argument('--horizon-log',action='store_true')
    parser.add_argument('--horizon-order',type=int)
    parser.add_argument('--reuse-source',type=Path)
    parser.add_argument('--extend-metric-source',type=Path,
                        help='Reuse lower ellmax source and compute only additional metric modes at matching radii')
    args=parser.parse_args()
    if min(args.radial_order,args.angular_order)<2:
        raise ValueError('Quadrature orders must be at least two')
    if args.background=='threshold-kerr':
        cloud=ThresholdCloud(alpha=args.alpha)
    else:
        from environment_schwarzschild_cloud import SchwarzschildCloud
        cloud=SchwarzschildCloud(alpha=args.alpha,freeze_decay=True)
    orbit=args.orbital_radius
    if orbit<=cloud.rp:
        raise ValueError('Orbit must be outside the horizon')
    metric=LorenzMetricMode(orbit,cloud.a,args.metric_m,args.metric_ellmax)
    scalar_m=cloud.m+args.metric_m
    scalar_ell=abs(scalar_m) if args.scalar_ell is None else args.scalar_ell
    if scalar_ell<abs(scalar_m):
        raise ValueError('Scalar ell must be >= |scalar m|')
    omega=cloud.omega+metric.omega
    green=RadialGreen(cloud.a,cloud.mu,omega,scalar_ell,scalar_m,
                      rmax=args.green_outer_radius,offset=args.green_horizon_offset,rtol=1e-11,
                      infinity_method=args.infinity_method)
    if args.infinity_method=='series' and green.series_last_term_relative>1e-3:
        raise ValueError('Unresolved inverse-r infinity boundary: use --infinity-method coulomb '
                         'and verify --green-outer-radius convergence before interpreting fluxes')
    inner=cloud.rp+args.source_inner_offset
    outer=args.source_outer_radius
    if not green.rmin<=inner<orbit<outer<=min(green.rmax,cloud.rmax):
        raise ValueError('Source cutoffs must enclose the orbit and lie inside both solved domains')
    panels=np.array([inner]+sorted({r for r in (3.,6.,12.,20.,40.,80.,160.,320.,orbit)
                                  if inner<r<outer})+[outer])
    node_panels,weight_panels=[],[]
    for index,(a,b) in enumerate(zip(panels[:-1],panels[1:])):
        order=(args.horizon_order or args.radial_order) if index==0 else args.radial_order
        if order<2:
            raise ValueError('Horizon quadrature order must be at least two')
        x,w=np.polynomial.legendre.leggauss(order)
        if index==0 and args.horizon_log:
            lo,hi=np.log(a-cloud.rp),np.log(b-cloud.rp)
            distance=np.exp((lo+hi)/2+(hi-lo)*x/2)
            node_panels.append(cloud.rp+distance)
            weight_panels.append((hi-lo)*w*distance/2)
        else:
            node_panels.append((a+b)/2+(b-a)*x/2)
            weight_panels.append((b-a)*w/2)
    radii,weights=np.concatenate(node_panels),np.concatenate(weight_panels)
    directory=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    suffix=''
    if args.alpha!=.3 or orbit!=20.:
        suffix=f'_alpha{args.alpha:g}_rp{orbit:g}'
    if args.background!='threshold-kerr':
        suffix+='_schwarzschild_frozen'
    if args.metric_m!=2 or scalar_ell!=3:
        suffix+=f'_mg{args.metric_m}_sl{scalar_ell}'
    if args.green_horizon_offset!=1e-4 or args.green_outer_radius!=1000.:
        suffix+=f'_gh{args.green_horizon_offset:g}_go{args.green_outer_radius:g}'
    if args.infinity_method!='series':
        suffix+='_'+args.infinity_method
    if args.source_inner_offset!=.05 or outer!=320.:
        suffix+=f'_inner{args.source_inner_offset:g}_outer{outer:g}'
    if args.horizon_log:
        suffix+='_log'
    if args.horizon_order is not None:
        suffix+=f'_h{args.horizon_order}'
    out=directory/f'forced_mode_nr{args.radial_order}_nt{args.angular_order}_L{args.metric_ellmax}{suffix}.json'
    metadata=dict(alpha=args.alpha,cloud_mass=1.,metric=metric.provenance,scalar_ell=scalar_ell,scalar_m=scalar_m,
                  omega=float(omega),radial_order=args.radial_order,angular_order=args.angular_order,
                  source_panels=panels.tolist(),green_outer_radius=args.green_outer_radius,
                  green_horizon_offset=args.green_horizon_offset)
    if args.infinity_method!='series':
        metadata['infinity_method']=args.infinity_method
    if args.background!='threshold-kerr':
        metadata['background']=dict(type=args.background,spectral_omega=encode(cloud.spectral_omega),
            normalization=cloud.normalization,
            approximation='Complex radial eigenfunction retained; temporal decay frozen; KG defect nonzero',
            source_operator='h^{ab} Hessian_ab as in paper; no claim of exact stationary balance')
    if args.horizon_log:
        metadata['horizon_log_first_panel']=True
    if args.horizon_order is not None:
        metadata['horizon_quadrature_order']=args.horizon_order
    reused={}
    extended={}
    additional_metric=None
    grid_keys={'radial_order','source_panels','green_outer_radius','green_horizon_offset',
               'horizon_log_first_panel','horizon_quadrature_order','infinity_method'}
    source_parameters=lambda p:{k:v for k,v in p.items() if k not in grid_keys}
    if args.reuse_source:
        cached=json.loads(args.reuse_source.read_text())
        if source_parameters(cached['parameters'])!=source_parameters(metadata):
            raise ValueError('Reusable samples have different source physics or angular truncation')
        reused={row['r']:row['source'] for row in cached['samples']}
    if args.extend_metric_source:
        cached=json.loads(args.extend_metric_source.read_text())
        lower=source_parameters(cached['parameters'])
        upper=source_parameters(metadata)
        lower=json.loads(json.dumps(lower))
        old_lmax=lower['metric']['ellmax']
        lower['metric']['ellmax']=args.metric_ellmax
        if lower!=upper or not abs(args.metric_m)<=old_lmax<args.metric_ellmax:
            raise ValueError('Extension requires identical source physics/angular quadrature and a lower complete metric ell sum')
        extended={row['r']:row['source'] for row in cached['samples']}
        additional_metric=LorenzMetricMode(orbit,cloud.a,args.metric_m,args.metric_ellmax,ellmin=old_lmax+1)
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
        if r in reused:
            value=reused[r]
        else:
            selected_metric=additional_metric if r in extended else metric
            frequency,J=project_source(cloud,[r],orbit,scalar_ell,scalar_m,selected_metric,ntheta=args.angular_order)
            if abs(frequency-omega)>1e-14:
                raise RuntimeError('Source and Green frequencies differ')
            value=encode(J[0]+(complex(*extended[r]) if r in extended else 0j))
        samples.append(dict(r=r,weight=float(weights[index]),source=value,reused=r in reused))
        if r in extended and r not in reused:
            samples[-1]['extended_from_metric_ellmax']=old_lmax
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
        cutoff_sequence.append(dict(outer_source_cutoff=float(cutoff),z_inf=encode(ci if green.propagating else 0j),
                                    up_coefficient=encode(ci),z_h=encode(ch)))
    for r in np.concatenate(([green.rmin],panels,[green.rmax])):
        left=radii<r
        cu=np.sum(weights[left]*ingoing[left]*J[left])/green.w0
        cv=np.sum(weights[~left]*outgoing[~left]*J[~left])/green.w0
        u,du=green.insol.sol(r)
        v,dv=green.upsol.sol(r)
        field_samples.append(dict(r=float(r),field=encode(v*cu+u*cv),derivative=encode(dv*cu+du*cv)))
    result.update(status='truncated_single_mode_not_converged',z_inf=encode(zi if green.propagating else 0j),
                  up_coefficient=encode(zi),propagating=bool(green.propagating),z_h=encode(zh),
                  radial_response_at_panel_boundaries=field_samples,
                  outer_source_cutoff_sequence=cutoff_sequence,
                  flux_scaling='per q^2*(cloud mass/M); unit Killing cloud mass, not alpha^-3 rescaled',
                  flux=mode_flux(omega,scalar_m,cloud.omega,cloud.m,cloud.mu,cloud.a,zi,zh),
                  wronskian_relative_spread=float(np.max(abs(green.wronskian(radii)/green.w0-1))),
                  infinity_boundary_audit=dict(method=green.infinity_method,
                      inverse_r_series_last_term_relative=green.series_last_term_relative,
                      status=('series_unresolved' if green.infinity_method=='series' and green.series_last_term_relative>1e-3
                              else 'outer_cutoff_convergence_required')),
                  missing_convergence=['radial quadrature','horizon source cutoff','outer source cutoff',
                                       'angular projection','metric ell truncation','remaining m modes and static completion'])
    save()
    print(json.dumps({k:result[k] for k in ('status','z_inf','z_h','flux','wronskian_relative_spread')}),flush=True)


if __name__=='__main__':
    main()
