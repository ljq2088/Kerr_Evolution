"""Bounded radial cross-backend audit, without modifying production selection."""
from __future__ import annotations
import argparse, hashlib, importlib.metadata, json, os, subprocess, sys, time
from pathlib import Path
import numpy as np
from pybhpt.radial import (RadialTeukolsky, renormalized_angular_momentum,
                           renormalized_angular_momentum_monodromy)
ROOT = Path(__file__).resolve().parents[1]
A = .8771530275949366
R0 = 20.
RADII = np.array([1 + np.sqrt(1 - A*A) + 1.e-3, 2., 10., 20., 40., 120., 316.])


def strict_json_values(value):
    """Keep failed numerical values explicit while emitting standard JSON."""
    if isinstance(value, dict):
        return {key: strict_json_values(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [strict_json_values(item) for item in value]
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return None
    return value


def enc(x):
    a = np.asarray(x, dtype=complex)
    return np.stack((a.real, a.imag), axis=-1).tolist()


def dec(x):
    a = np.asarray(x, dtype=float)
    return a[..., 0] + 1j*a[..., 1]


def child(case):
    start = time.monotonic()
    s, ell, m = case['spin'], case['ell'], case['m']
    om = m/(R0**1.5 + A)
    radial = RadialTeukolsky(s, ell, m, A, om, RADII)
    if case['method'] == 'NU':
        nu = renormalized_angular_momentum(s, ell, m, A, om)
        mono = renormalized_angular_momentum_monodromy(s, ell, m, A, om, radial.eigenvalue)
        out = dict(case, status='ok', omega=om, eigenvalue=radial.eigenvalue,
                   selected_nu=enc(nu), monodromy_nu=enc(mono),
                   selected_nu_finite=bool(np.isfinite(nu)),
                   selected_nu_abs_difference_from_monodromy=float(abs(nu-mono)),
                   elapsed_seconds=time.monotonic()-start)
        print('RESULT_JSON='+json.dumps(out, allow_nan=True), flush=True)
        return
    radial.solve(method=case['method'], rtol=case.get('rtol'))
    val = np.array([radial(bc, d) for bc in ('In', 'Up') for d in (0, 1)]).reshape(2, 2, -1)
    finite = bool(np.isfinite(val).all())
    out = dict(case, status='ok' if finite else 'nonfinite', omega=om,
               elapsed_seconds=time.monotonic()-start, eigenvalue=radial.eigenvalue,
               values=enc(val))
    if np.isfinite(val[:,:,3]).all():
        delta = RADII**2 - 2*RADII + A*A
        w = delta**(s+1) * (val[0,0]*val[1,1] - val[1,0]*val[0,1])
        wref = w[3]
        kernel = np.where(RADII < R0, val[1,0,3]*val[0,0], val[0,0,3]*val[1,0])/wref
        dkernel = np.where(RADII < R0, val[1,0,3]*val[0,1], val[0,0,3]*val[1,1])/wref
        out.update(wronskian=enc(w), kernel=enc(kernel), kernel_derivative=enc(dkernel),
                   wronskian_max_relative_drift=float(np.max(abs(w/wref - 1))),
                   boundary_points=[float(radial.boundarypoint(bc)) for bc in ('In','Up')],
                   boundary_values=enc([radial.boundarysolution(bc) for bc in ('In','Up')]),
                   boundary_derivatives=enc([radial.boundaryderivative(bc) for bc in ('In','Up')]))
    print('RESULT_JSON='+json.dumps(out, allow_nan=True), flush=True)


def compare(row, base):
    if 'values' not in row or base['status'] != 'ok': return
    v, b = dec(row['values']), dec(base['values'])
    row['comparison_reference'] = dict(method=base['method'], rtol=base.get('rtol'))
    physical = np.where(RADII[None,:] <= R0, v[0], v[1])
    physical_base = np.where(RADII[None,:] <= R0, b[0], b[1])
    row['physical_branch_finite'] = bool(np.isfinite(physical).all())
    row['physical_branch_raw_relative_by_derivative_and_radius'] = abs(physical/physical_base-1).tolist()
    row['physical_branch_raw_max_relative_change'] = float(np.max(abs(physical/physical_base-1)))
    row['physical_branch_raw_max_relative_change_r_ge_2'] = float(np.max(abs(physical[:,1:]/physical_base[:,1:]-1)))
    row['raw_max_relative_change'] = float(np.max(abs(v/b - 1)))
    row['raw_max_relative_by_solution_derivative'] = np.max(abs(v/b-1), axis=-1).tolist()
    row['normalization_ratio_at_r20'] = enc(v[:,0,3]/b[:,0,3])
    row['constant_rescaled_shape_max_relative_change'] = float(np.max(abs(v / (v[:,0,3]/b[:,0,3])[:,None,None] / b - 1)))
    row['logarithmic_derivative_max_relative_change'] = float(np.max(abs((v[:,1]/v[:,0])/(b[:,1]/b[:,0])-1)))
    if 'kernel' in row:
        row['green_kernel_max_relative_change'] = float(np.max(abs(dec(row['kernel'])/dec(base['kernel'])-1)))
        row['green_kernel_derivative_max_relative_change'] = float(np.max(abs(dec(row['kernel_derivative'])/dec(base['kernel_derivative'])-1)))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--child'); p.add_argument('--timeout', type=float, default=25.)
    p.add_argument('--group', choices=['dominant','frequency','highell','precision','nu'], default='dominant')
    p.add_argument('--output', default=None)
    args = p.parse_args()
    if args.child: child(json.loads(args.child)); return
    if args.group == 'dominant':
        modes = [(s,max(1,abs(s)),1) for s in (0,-1,1,-2,2)]
    elif args.group == 'frequency':
        modes = [(s,max(2,m),m) for m in range(2,7) for s in (-2,2)]
    elif args.group == 'highell':
        modes = [(s,ell,1) for ell in (6,18) for s in (0,-1,1,-2,2)]
    elif args.group == 'precision':
        modes = [(-2,2,1),(-2,3,3)]
    else:
        modes = [(0,l,1) for l in range(1,7)] + [(s,l,1) for s in (-1,1) for l in (1,6,18)] + [(-2,2,1),(2,2,1)]
    cases = []
    for s,ell,m in modes:
        choices = [('AUTO',None),('HBL',1.e-15),('MST',None)] + ([('GSN',1.e-14)] if abs(s)==2 else [])
        if args.group == 'precision':
            choices = [('AUTO',None),('HBL',1.e-15),('TEUK',1.e-14)] + [('GSN',rtol) for rtol in (1.e-11,1.e-12,1.e-13,1.e-14,1.e-15)]
        elif args.group == 'nu':
            choices = [('NU',None)]
        for method,rtol in choices:
            cases.append(dict(spin=s,ell=ell,m=m,method=method,rtol=rtol))
    src = ROOT/'outputs/lorenz_reference/pybhpt_radial_source/v1.0.0/radialsolver.cpp'
    import cybhpt_full
    out = dict(status='bounded_radial_comparison_not_full_flux_validation', a=A, rp=R0,
      pybhpt_version=importlib.metadata.version('pybhpt'), radii=RADII.tolist(),
      binary_sha256=hashlib.sha256(Path(cybhpt_full.__file__).read_bytes()).hexdigest(),
      report_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      reference_source_commit='9fe9c57e2d1c92d944ba70ba3c1b81b665e9d274',
      reference_source_url='https://github.com/znasipak/pybhpt/blob/9fe9c57e2d1c92d944ba70ba3c1b81b665e9d274/cpp/src/radialsolver.cpp',
      reference_source_sha256=hashlib.sha256(src.read_bytes()).hexdigest() if src.exists() else None,
      nonfinite_json_encoding='null denotes a nonfinite or undefined numerical value; inspect finite/status flags',
      normalization='Common pybhpt In/Up asymptotic normalization; raw and invariant Green products compared.',
      default_ode_relative_tolerance=1.e-14,
      actual_algorithm_from_versioned_source={
        'AUTO':'ASYMP only if abs(omega*r_min)>100; otherwise HBL with stable spin flip. MST block commented out. No auto GSN.',
        'MST':'Direct analytic MST series; not guaranteed to succeed for every mode/radius.',
        'GSN':'Only abs(spin)=2. Positive spin uses spin flip.',
        'HBL':'Hyperboloidally transformed ODE. In s>0 and Up s<0 use spin flip.'},
      limitations=['No amplitude/source reconstruction or scalar flux calculated here.',
        'Wronskian constancy does not check asymptotic amplitude normalization.',
        'Cross-method agreement is evidence, not an exact true-error estimate.',
        'Boundary generation is shared by HBL/GSN; MST is independent of ODE propagation.',
        'GSN is not called for unsupported spin 0 or 1.',
        'Selection inferred from exact version tag; wrapper bytes match installed version. Binary branch not instrumented.'], rows=[])
    dest = Path(args.output) if args.output else ROOT/f'docs/environment_reproduction/radial_backend_audit_{args.group}.json'
    for case in cases:
        start = time.monotonic()
        try:
            run = subprocess.run([sys.executable, str(Path(__file__).resolve()), '--child', json.dumps(case)],
                capture_output=True, text=True, timeout=args.timeout, env={**os.environ,'OPENBLAS_NUM_THREADS':'1'})
            marker = next((x for x in run.stdout.splitlines() if x.startswith('RESULT_JSON=')),None)
            row = json.loads(marker[len('RESULT_JSON='):]) if marker else dict(case,status='failed',returncode=run.returncode)
            row['diagnostic_output'] = '\n'.join(x for x in run.stdout.splitlines() if not x.startswith('RESULT_JSON='))[-4000:]
            row['stderr'] = run.stderr[-2000:]
        except subprocess.TimeoutExpired:
            row = dict(case,status='timeout',elapsed_seconds=time.monotonic()-start)
        base = next((x for x in out['rows'] if all(x[k]==case[k] for k in ('spin','ell','m')) and x['method']=='AUTO'),None)
        if base: compare(row,base)
        out['rows'].append(row)
        tmp = dest.with_suffix('.tmp'); tmp.write_text(json.dumps(strict_json_values(out),indent=2,allow_nan=False)+'\n');tmp.replace(dest)
        print(json.dumps({k:v for k,v in row.items() if k not in ('values','wronskian','kernel','kernel_derivative','boundary_values','boundary_derivatives','diagnostic_output','stderr','physical_branch_raw_relative_by_derivative_and_radius')}),flush=True)

if __name__=='__main__': main()
