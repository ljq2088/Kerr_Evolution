"""Audit plotted radial interpolation against the saved Gauss Green solution."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import json
from pathlib import Path
import numpy as np
from environment_response import SampledResponse
from environment_source import angular_mode


def audit(path):
    data=json.loads(path.read_text());p=data['parameters']
    response=SampledResponse.from_report(data)
    old=np.array([complex(*r['field']) for r in data['radial_response_at_panel_boundaries']])
    radii=np.array([r['r'] for r in data['radial_response_at_panel_boundaries']])
    new=response.evaluate(radii)[0]
    amp={}
    for key,value in [('z_h',response.horizon_coefficient),('up_coefficient',response.up_coefficient)]:
        before=complex(*data[key])
        amp[key]=dict(relative_difference=float(abs(value-before)/abs(before)) if before else None,
            saved=[before.real,before.imag],continuous=[value.real,value.imag])
    probe=np.geomspace(max(response.green.rmin,response.green.rp+.05),318.5,320)
    field=response.evaluate(probe)[0]
    S=angular_mode(np.array([np.pi/2]),p['scalar_ell'],p['scalar_m'],
        p['metric']['a']**2*(p['omega']**2-p['alpha']**2))[0][0]
    scaled=abs(field*S*p['alpha']**-3);index=int(np.argmax(scaled))
    if not np.isfinite(new).all() or not np.isfinite(scaled).all():raise ValueError('Nonfinite response')
    return dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        ell=p['scalar_ell'],m=p['scalar_m'],amplitudes=amp,
        panel_field_difference_over_saved_max=float(np.max(abs(new-old))/np.max(abs(old))) if np.max(abs(old)) else None,
        equatorial_mode_maximum=float(scaled[index]),maximum_r=float(probe[index]))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('coverage',type=Path);parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.workers<1:raise ValueError('Workers must be positive')
    data=json.loads(args.coverage.read_text());modes=data['field']['modes']
    if any('flux' not in row for row in modes):raise ValueError('Require complete field inventory')
    paths=[args.coverage.parent/row['file'] for row in modes]
    result=dict(status='radial_consistency_in_progress',coverage=str(args.coverage),
        coverage_sha256=hashlib.sha256(args.coverage.read_bytes()).hexdigest(),expected=len(paths),rows=[])
    def save():
        temporary=args.output.with_suffix('.tmp');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(args.output)
    save()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(audit,p) for p in paths]
        for future in as_completed(futures):
            row=future.result();result['rows'].append(row);save()
            print(row['ell'],row['m'],row['panel_field_difference_over_saved_max'],row['equatorial_mode_maximum'],flush=True)
    result['rows'].sort(key=lambda row:(row['ell'],row['m']))
    result['status']='completed_radial_consistency_not_convergence'
    result['limitation']='Same saved source nodes; interpolated-source integration versus original Gauss integration. Neither alone proves accuracy against the continuum source or paper.'
    save()


if __name__=='__main__':main()
