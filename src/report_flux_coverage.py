"""Inventory a common finite-resolution flux set without summing missing modes."""
import argparse
import json
from datetime import datetime,timezone
from pathlib import Path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--metric-ellmax',type=int,default=6)
    parser.add_argument('--orbital-radius',type=float,default=20.)
    parser.add_argument('--angular-order',type=int,default=10)
    parser.add_argument('--infinity-ellmax',type=int,default=6)
    parser.add_argument('--horizon-ellmax',type=int,default=5)
    parser.add_argument('--field-ellmax',type=int,default=5)
    parser.add_argument('--green-outer-by-m',nargs=2,action='append',default=[],
                        metavar=('SCALAR_M','RADIUS'),
                        help='Explicit Green outer radius for one scalar m sector; default 1000M')
    parser.add_argument('--infinity-method-by-m',nargs=2,action='append',default=[],
                        metavar=('SCALAR_M','METHOD'),
                        help='Explicit series or coulomb boundary construction for one scalar m sector')
    args=parser.parse_args()
    from environment_cloud import cloud_211
    cloud=cloud_211(outer_efolds=45.,horizon_offset=1e-6,rtol=2e-11,mu=.3)
    a=cloud['a_over_M'];omega_c=cloud['M_omega_c']
    if args.orbital_radius<=cloud['r_plus_over_M']:
        raise ValueError('Orbit must lie outside the horizon')
    omega_p=1/(args.orbital_radius**1.5+a)
    boundaries={}
    for m,radius in args.green_outer_by_m:
        m=int(m);radius=float(radius)
        if m in boundaries or radius<=320:raise ValueError('Invalid or duplicate boundary override')
        boundaries[m]=radius
    methods={}
    for m,method in args.infinity_method_by_m:
        m=int(m)
        if m in methods or method not in ('series','coulomb'):
            raise ValueError('Invalid or duplicate infinity method override')
        methods[m]=method
    folder=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    found={}
    for path in folder.glob('forced_mode_*.json'):
        data=json.loads(path.read_text());p=data.get('parameters',{});g=p.get('metric',{})
        # Alternative angular diagnostics are separate convergence experiments.
        if g.get('angular_backend') is not None:
            continue
        panels=p.get('source_panels',[])
        if (p.get('alpha')!=.3 or g.get('orbital_radius')!=args.orbital_radius or p.get('background') is not None
            or abs(g.get('a',0)-a)>1e-12
            or g.get('ellmax')!=args.metric_ellmax or p.get('radial_order')!=8
            or p.get('angular_order')!=args.angular_order or p.get('horizon_quadrature_order')!=32
            or not p.get('horizon_log_first_panel') or not panels
            or abs(panels[0]-(1+(1-g['a']**2)**.5)-.0005)>1e-12
            or panels[-1]!=320 or p.get('green_outer_radius')!=boundaries.get(p.get('scalar_m'),1000)
            or p.get('cloud_mass')!=1
            or p.get('green_horizon_offset')!=.0001 or p.get('infinity_method','series')!=methods.get(p.get('scalar_m'),'series')):
            continue
        key=(p['scalar_ell'],p['scalar_m'])
        if key in found:raise ValueError(f'Ambiguous duplicate mode {key}')
        found[key]=dict(file=path.name,status=data['status'],samples=len(data['samples']),
                       green_outer_radius=p['green_outer_radius'],infinity_method=p.get('infinity_method','series'))
        if (data['status']=='truncated_single_mode_not_converged' and 'flux' in data
            and data.get('flux_validity')!='historical_unreliable_boundary_result'):
            found[key]['flux']=data['flux']
    # Equatorial reflection symmetry requires ell+m even for the |211> cloud.
    needed_infinity=[(ell,m) for ell in range(args.infinity_ellmax+1) for m in range(-ell,ell+1)
                     if (ell+m)%2==0 and m!=1 and (omega_c+(m-1)*omega_p)**2>.3**2]
    needed_horizon=[(ell,m) for ell in range(args.horizon_ellmax+1) for m in range(-ell,ell+1)
                    if (ell+m)%2==0 and m!=1]
    needed_field=[(ell,m) for ell in range(2,args.field_ellmax+1) for m in range(-ell,ell+1)
                  if (ell+m)%2==0]
    field_rows=[dict(ell=ell,m=m,**found.get((ell,m),dict(status='missing'))) for ell,m in needed_field]
    def boundary(name,needed):
        rows=[dict(ell=ell,m=m,**found.get((ell,m),dict(status='missing'))) for ell,m in needed]
        complete=all('flux' in row for row in rows)
        return dict(required_count=len(rows),computed_count=sum('flux' in row for row in rows),
            modes=rows,finite_resolution_total=(sum(row['flux'][name]['orbital_energy'] for row in rows) if complete else None),
            converged=False)
    report=dict(observed_utc=datetime.now(timezone.utc).isoformat(),
        status='coverage_only_not_paper_reproduction',
        parameters=dict(alpha=.3,rp=args.orbital_radius,metric_ellmax=args.metric_ellmax,nr=8,nt=args.angular_order,horizon_order=32,
                        inner_offset=.0005,outer_source_cutoff=320,green_outer_by_scalar_m=boundaries,infinity_method_by_scalar_m=methods),
        infinity=boundary('infinity',needed_infinity),horizon=boundary('horizon',needed_horizon),
        field=dict(required_count=len(field_rows),computed_count=sum('flux' in row for row in field_rows),
                   modes=field_rows,converged=False),
        excluded='m=1 is static forcing: orbital-effective energy flux vanishes, but its field is required for wakes. Nonpropagating infinity modes carry zero stationary infinity flux.',
        limits='No total is reported with missing modes. Even a complete finite-resolution total needs L_g=18 and all other convergence checks.')
    suffix=('' if args.angular_order==10 and args.infinity_ellmax==6 and args.horizon_ellmax==5
            and args.field_ellmax==5 and not boundaries and not methods else
            f'_nt{args.angular_order}_i{args.infinity_ellmax}_h{args.horizon_ellmax}_f{args.field_ellmax}')
    if boundaries:suffix+='_'+ '_'.join(f'm{m}go{radius:g}' for m,radius in sorted(boundaries.items()))
    if methods:suffix+='_'+ '_'.join(f'm{m}{method}' for m,method in sorted(methods.items()))
    if args.orbital_radius!=20.:suffix+=f'_rp{args.orbital_radius:g}'
    (folder/f'flux_coverage_L{args.metric_ellmax}{suffix}.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:{a:report[k][a] for a in ('required_count','computed_count','finite_resolution_total')} for k in ('infinity','horizon')}))


if __name__=='__main__':main()
