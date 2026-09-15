"""Prepare local source-jump references and compare completed original Sam output.

Sam's m=1 closure uses l-l- value/derivative. The production file is not
modified: an isolated function copy changes only that row selection and
returns its assembled matrix for a cheap same-matrix mixed-row control.
"""
import argparse
import hashlib
import inspect
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def encode(x):
    a = np.asarray(x, complex)
    return np.stack([a.real, a.imag], axis=-1).tolist()


def decode(x):
    a = np.asarray(x, float)
    return a[..., 0] + 1j*a[..., 1]


def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(data, indent=2)+'\n')
    tmp.replace(path)


def prepare_matching(args):
    from lorenz_jet import Jet
    from paper_precise_angular import install_matching_angular
    from paper_analytic_jumps import curvature_jumps
    import paper_all_component_matching as matching
    import mpmath as mp
    Jet.coefficient_dtype = np.clongdouble
    install_matching_angular()
    matching.curvature_jumps = curvature_jumps
    source = inspect.getsource(matching.solve)
    old = 'if m==1:mask[:,0,2]=False;mask[:,0,4]=True'
    new = 'if m==1:mask[:,0,2]=False;mask[:,0,1]=True'
    if source.count(old) != 1 or source.count('return dict(r0=') != 1:
        raise RuntimeError('Matching implementation changed; audit the isolated transformation')
    amended = source.replace(old, new).replace('return dict(r0=', 'return dict(operator_matrix=A,known_projected=known_proj,r0=')
    namespace = dict(vars(matching))
    exec(compile(amended, '<isolated_original_sam_row_selection>', 'exec'), namespace)
    result = namespace['solve'](args.r0, args.a, args.m, args.ellmax, args.quadrature, extended_solve=True)
    matrix, known = result.pop('operator_matrix'), result.pop('known_projected')
    mask = result['selected_mask'].copy()
    if args.m == 1:
        mask[:, 0, 1] = False
        mask[:, 0, 4] = True
    B = matrix[mask]
    b = (result['target']-known)[mask]
    rs = np.linalg.norm(B, axis=1)
    B, b = B/rs[:,None], b/rs
    cs = np.linalg.norm(B, axis=0)
    B = B/cs[None,:]
    def convert(z):
        fmt = lambda v: np.format_float_scientific(np.longdouble(v), precision=35, unique=False)
        return mp.mpc(fmt(np.real(z)), fmt(np.imag(z)))
    with mp.workdps(50):
        x = mp.lu_solve(mp.matrix([[convert(z) for z in row] for row in B]), mp.matrix([convert(z) for z in b]))
        x = np.array([np.clongdouble(np.longdouble(str(z.real)))+1j*np.longdouble(str(z.imag)) for z in x])/cs
    control_residual = np.einsum('djfc,c->djf', matrix, x)+known-result['target']
    result['mixed_row_control_solution'] = x
    result['mixed_row_control_residual'] = control_residual
    result['matching_row_choice'] = 'Sam m=1: l+l+ and m+m+ except ell1, then l-l- value/derivative'
    result['same_matrix_control'] = 'Current diagnostic m=1: replace ell1 l-l- with rho l+m+ value/derivative'
    result['assembly_function_sha256'] = hashlib.sha256(source.encode()).hexdigest()
    result['isolated_transformed_function_sha256'] = hashlib.sha256(amended.encode()).hexdigest()
    for key, value in list(result.items()):
        if isinstance(value, np.ndarray):
            result[key] = value.tolist() if value.dtype == bool else encode(value)
    result['status'] = 'local_sam_row_matching_prepared'
    result['trusted_degrees'] = [l for l in range(args.m, args.ellmax+1) if args.ellmax-l > 6]
    result['limitations'] = ['Original author data not used in this local assembly',
        '50-digit linear solve over extended precision inputs, not 50-digit physical accuracy',
        'Highest seven degrees retained as edge diagnostics only; no total flux comparison']
    from source_provenance import local_dependency_hashes
    result['implementation_sha256'] = local_dependency_hashes(ROOT/'src', ['paper_all_component_matching', 'paper_precise_angular', 'paper_analytic_jumps'])
    save(args.output, result)


def prepare_production(args):
    from environment_angular_diagnostic import install_dense_angular_diagnostic
    from lorenz_jet import Jet
    from report_paper_dipole_closure import old_spin1_jumps
    from report_paper_kappa_tables import calculate
    from source_provenance import local_dependency_hashes
    install_dense_angular_diagnostic()
    Jet.coefficient_dtype = np.complex128
    result = dict(status='local_production_mapped_jumps', a=args.a, r0=args.r0, m=args.m, L=args.ellmax,
        unknown_order=['Pm1j0','Pm1j1','kappaj0','kappaj1'], rows=[],
        implementation_sha256=local_dependency_hashes(ROOT/'src', ['report_paper_dipole_closure','report_paper_kappa_tables']),
        limitations=['Spin1 jumps extracted from production full-current gauge vector with quadratic orbit extrapolation and recorded basis-fit residual',
          'Kappa23 mapping includes angular source-projector derivative and chi; not direct kappa24',
          'No original author data consumed and no environmental flux evaluated'])
    for ell in range(args.m, args.ellmax+1):
        vector, residual = old_spin1_jumps(args.r0, args.a, args.m, ell)
        k0, k1, correction = calculate(args.a, args.r0, args.m, ell)
        result['rows'].append(dict(ell=ell,jumps=encode(np.r_[vector,k0,k1]),
            spin1_basis_fit_relative_residual=residual,kappa_angular_contact_correction=float(correction)))
        save(args.output,result)
        print('mapped production jumps',ell,'fit residual',residual,flush=True)


def compare(args):
    # No polling: consume only a complete valid file supplied explicitly.
    author = json.loads(args.author.read_text())
    local = json.loads(args.local.read_text())
    for key in ('a','r0','m'):
        if not np.isclose(author[key],local[key],rtol=0,atol=1e-14):
            raise ValueError('Author/local physical parameter mismatch: '+key)
    if author['unknown_order'] != ['Pm1j0','Pm1j1','kappaj0','kappaj1']:
        raise ValueError('Unexpected author jump order')
    ref = decode(author['jumps'])
    if ref.shape != (author['lmax'],4) or author['m'] != 1:
        raise ValueError('This author driver exports ell=1..L, m=1 only')
    if local['status']=='local_production_mapped_jumps':
        values={row['ell']:decode(row['jumps']) for row in local['rows']}
        mixed={}
    else:
        labels=local['labels'];vals=decode(local['solution']);ctrl=decode(local['mixed_row_control_solution'])
        values={};mixed={}
        for i,(ell,kind,datum) in enumerate(labels):
            k={'spin1':0,'kappa':2}[kind]+datum
            values.setdefault(ell,np.zeros(4,complex))[k]=vals[i]
            mixed.setdefault(ell,np.zeros(4,complex))[k]=ctrl[i]
    rows=[]
    for ell,value in sorted(values.items()):
        if ell>author['lmax']:continue
        reference=ref[ell-1]
        difference=value-reference
        scale=max(float(np.max(abs(reference))),1e-30)
        row=dict(ell=ell,trusted_after_six_mode_guard=author['lmax']-ell>6,
            author_jumps=encode(reference),local_jumps=encode(value),difference=encode(difference),
            absolute_max=float(np.max(abs(difference))),row_scale=scale,relative_row_max=float(np.max(abs(difference))/scale))
        if ell in mixed:row['mixed_row_control_difference']=encode(mixed[ell]-reference)
        rows.append(row)
    save(args.output,dict(status='completed_original_sam_jump_comparison',
        author_sha256=hashlib.sha256(args.author.read_bytes()).hexdigest(),
        local_sha256=hashlib.sha256(args.local.read_bytes()).hexdigest(),
        author_lmax=author['lmax'],local_lmax=local['L'],rows=rows,
        limitations=['Relative scale is each four-entry jump row; symmetry-zero entries are not assigned misleading relative error',
            'Ell4..10 of an L10 calculation are edge diagnostics, not validated interior modes',
            'Comparison of local jump inputs; not a final environmental flux test']))


def main():
    p=argparse.ArgumentParser()
    p.add_argument('mode',choices=['matching','production','compare'])
    p.add_argument('--a',type=float,default=.6);p.add_argument('--r0',type=float,default=8.)
    p.add_argument('--m',type=int,default=1);p.add_argument('--ellmax',type=int,default=10)
    p.add_argument('--quadrature',type=int,default=32)
    p.add_argument('--author',type=Path);p.add_argument('--local',type=Path)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    {'matching':prepare_matching,'production':prepare_production,'compare':compare}[args.mode](args)


if __name__=='__main__':main()
