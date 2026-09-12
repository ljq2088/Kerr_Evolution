"""Plan and optionally compute missing nonstatic channels from a flux inventory.

Groups opposite metric m branches only when their Green settings agree.
Static channels require their separate matching problem and are never skipped.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


COMPLETE = 'truncated_single_mode_not_converged'


def build_plan(coverage_path, targets=('field','infinity','horizon'), workers=2):
    if workers < 1:
        raise ValueError('Workers must be positive')
    data = json.loads(coverage_path.read_text())
    if data['status'] != 'coverage_only_not_paper_reproduction':
        raise ValueError('Expected a coverage inventory')
    p = data['parameters']
    needed = set()
    verified = {}
    for target in targets:
        if target not in ('field','infinity','horizon'):
            raise ValueError('Unknown target')
        for row in data[target]['modes']:
            mode = (row['ell'], row['m'])
            if mode in verified:
                continue
            complete = False
            if row.get('status') == COMPLETE and 'flux' in row and 'file' in row:
                path = coverage_path.parent / row['file']
                if path.is_file():
                    actual = json.loads(path.read_text())
                    q = actual['parameters']
                    if (q['scalar_ell'], q['scalar_m']) != mode:
                        raise ValueError('Inventory mode differs from response')
                    if (q['alpha'],q['metric']['orbital_radius'],q['metric']['a']) != (p['alpha'],p['rp'],p['a']):
                        raise ValueError('Inventory physics differs from response')
                    expected = (p['metric_ellmax'],p['nr'],p['nt'],p['horizon_order'],
                                p['green_outer_by_scalar_m'].get(str(mode[1]),1000.),
                                p['infinity_method_by_scalar_m'].get(str(mode[1]),'series'))
                    observed = (q['metric']['ellmax'],q['radial_order'],q['angular_order'],
                                q['horizon_quadrature_order'],q['green_outer_radius'],
                                q.get('infinity_method','series'))
                    if observed != expected:
                        raise ValueError('Inventory discretization differs from response')
                    complete = (actual.get('status') == COMPLETE and actual.get('flux') == row['flux']
                                and actual.get('flux_validity') != 'historical_unreliable_boundary_result')
                    if complete:
                        verified[mode] = dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            if not complete:
                needed.add(mode)
    needed.difference_update(verified)
    grouped = {}
    static = []
    for ell,m in sorted(needed):
        mg = m-1
        if mg == 0:
            static.append([ell,m])
            continue
        outer = p['green_outer_by_scalar_m'].get(str(m),1000.)
        method = p['infinity_method_by_scalar_m'].get(str(m),'series')
        group = grouped.setdefault((abs(mg),outer,method),dict(direct=[],conjugate=[]))
        group['direct' if mg > 0 else 'conjugate'].append(ell)
    jobs = []
    for (mg,outer,method),group in sorted(grouped.items()):
        args = ['--alpha',str(p['alpha']),'--background',p['background'],
                '--orbital-radius',str(p['rp']),'--metric-m',str(mg),
                '--metric-ellmax',str(p['metric_ellmax']),'--radial-order',str(p['nr']),
                '--angular-order',str(p['nt']),'--horizon-log',
                '--horizon-order',str(p['horizon_order']),
                '--source-inner-offset',str(p['inner_offset']),
                '--source-outer-radius',str(p['outer_source_cutoff']),
                '--green-horizon-offset','0.0001','--green-outer-radius',str(outer),
                '--infinity-method',method,'--workers',str(workers)]
        for key,flag in [('direct','--scalar-ells'),('conjugate','--conjugate-ells')]:
            if group[key]:args += [flag,*map(str,sorted(group[key]))]
        jobs.append(dict(metric_m=mg,scalar_ells=sorted(group['direct']),
                         conjugate_ells=sorted(group['conjugate']),arguments=args))
    return dict(status='missing_mode_plan_not_convergence',
        inventory=dict(file=str(coverage_path),sha256=hashlib.sha256(coverage_path.read_bytes()).hexdigest()),
        targets=list(targets),verified_existing_count=len(verified),
        missing_modes=[list(mode) for mode in sorted(needed)],
        unresolved_static_modes=static,jobs=jobs,
        limitation='Completes the requested finite inventory only. Static matching and convergence require separate checks.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('coverage',type=Path)
    parser.add_argument('--targets',nargs='+',choices=('field','infinity','horizon'),default=['field','infinity','horizon'])
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--execute',action='store_true')
    args = parser.parse_args()
    result = build_plan(args.coverage,args.targets,args.workers)
    def save():
        args.output.parent.mkdir(parents=True,exist_ok=True)
        tmp=args.output.with_suffix('.tmp');tmp.write_text(json.dumps(result,indent=2)+'\n');tmp.replace(args.output)
    save()
    if args.execute:
        if result['unresolved_static_modes']:
            raise ValueError('Static modes require explicit matching data; no jobs started')
        result['status']='missing_mode_execution_in_progress';save()
        for job in result['jobs']:
            command=[sys.executable,str(Path(__file__).with_name('report_environment_mode_batch.py')),*job['arguments']]
            outcome=subprocess.run(command,check=False)
            job['returncode']=outcome.returncode
            if outcome.returncode:
                result['status']='missing_mode_execution_failed';save()
                raise RuntimeError('Batch failed; completed channel outputs preserved')
            save()
        result['status']='batches_completed_refresh_inventory_required';save()
    print(json.dumps({k:result[k] for k in ('status','missing_modes','unresolved_static_modes')}),flush=True)


if __name__ == '__main__':
    main()
