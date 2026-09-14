"""Compare mapped scalar-potential jumps with three published 2023 tables.

The convention map is fixed before comparison; no amplitude is fitted.
"""
import json,re,hashlib
from pathlib import Path
import numpy as np
from environment_angular_diagnostic import install_dense_angular_diagnostic
from environment_angular_variation import angular_mode_mass_derivative
from environment_source import kerr_metric
from lorenz_chi import chi_amplitudes
from lorenz_metric import _homogeneous_radial_data


def published_tables():
    root=Path(__file__).resolve().parents[1]
    source=root/'outputs/lorenz_reference/LorenzGaugeKerrCirc.tex'
    tex=source.read_text();tables=[]
    for label,a in (('tbl:jumps-a0',0.),('tbl:jumps-a06',.6),('tbl:jumps-a099',.99)):
        end=tex.index('\\label{'+label+'}')
        start=tex.rfind('\\begin{table}',0,end)
        block=tex[start:end];data=[]
        for line in block.splitlines():
            if not re.match(r'^\s*\d+\s*&',line):continue
            cells=[s.strip().replace('\\','').strip() for s in line.split('&')]
            if len(cells)!=5:continue
            data.append(dict(ell=int(cells[0]),kappa0_text=cells[3],kappa1_text=cells[4],
                        kappa0=float(cells[3]),kappa1=float(cells[4])))
        if len(data)!=9:raise ValueError('Unexpected published table shape')
        tables.append(dict(a=a,r0=6.,m=2,label=label,rows=data))
    return tables,hashlib.sha256(source.read_bytes()).hexdigest()


def calculate(a,r0,m,ell):
    omega=m/(r0**1.5+a);op=omega/m;delta=r0*r0-2*r0+a*a
    g=kerr_metric(r0,np.pi/2,a)
    ut=1/np.sqrt(-g[0,0]-2*op*g[0,3]-op*op*g[3,3])
    _,_,_,ds,_,_=angular_mode_mass_derivative(np.pi/2,ell,m,a,omega,0.)
    up,inn=chi_amplitudes(r0,a,ell,m)
    _,ru,du=_homogeneous_radial_data(0,ell,m,a,omega,r0,'Up')
    _,ri,di=_homogeneous_radial_data(0,ell,m,a,omega,r0,'In')
    j0=-1j*(up*ru-inn*ri)/(2*omega)
    correction=-8*np.pi*ds/(ut*delta)
    j1=-1j*(up*du-inn*di)/(2*omega)+correction
    return j0,j1,correction


def main():
    install_dense_angular_diagnostic()
    tables,digest=published_tables();rows=[]
    for table in tables:
      for ref in table['rows']:
        v0,v1,correction=calculate(table['a'],table['r0'],table['m'],ref['ell'])
        row=dict(a=table['a'],r0=6.,m=2,ell=ref['ell'],reference_label=table['label'],
          reference_kappa0=ref['kappa0'],reference_kappa1=ref['kappa1'],
          calculated_kappa0=[v0.real,v0.imag],calculated_kappa1=[v1.real,v1.imag],
          angular_source_derivative_correction=float(correction),
          absolute_error0=float(abs(v0-ref['kappa0'])),absolute_error1=float(abs(v1-ref['kappa1'])))
        # Roundoff allowance is derived from each printed nonzero decimal,
        # with a separate absolute numerical floor for symmetry-forbidden zeros.
        passed=[]
        for k,value in ((0,v0),(1,v1)):
            literal=ref[f'kappa{k}_text'];reference=ref[f'kappa{k}']
            tolerance=(.5*10**(-len(literal.split('.')[1])) if '.' in literal else 0.)+1e-7
            row[f'comparison_tolerance{k}']=tolerance
            passed.append(abs(value-reference)<=tolerance)
        row['within_printed_precision_plus_1e_minus7_absolute']=bool(all(passed))
        rows.append(row);print(json.dumps(row),flush=True)
    result=dict(status='mapped_scalar_jumps_compared_to_published_tables',
      source_url='https://arxiv.org/html/2306.16459v3',source_sha256=digest,
      map='kappa23=i*(kappa24-chi_compact)/(2*omega)',
      derivative_map='J1(kappa_jj23)=-i*J1(chi_j)/(2*omega)-8*pi*d_mass_squared(S_j(pi/2))/(ut*Delta0)',
      all_rows_within_printed_precision_plus_1e_minus7_absolute=all(row['within_printed_precision_plus_1e_minus7_absolute'] for row in rows),
      maximum_nonzero_relative_error=max(row[f'absolute_error{k}']/abs(row[f'reference_kappa{k}']) for row in rows for k in (0,1) if row[f'reference_kappa{k}']!=0),
      maximum_forbidden_absolute_value=max(row[f'absolute_error{k}'] for row in rows for k in (0,1) if row[f'reference_kappa{k}']==0),
      rows=rows,limitations=['Published values have about 7-8 significant digits',
        'Comparison covers r0=6,m=2,a=0,0.6,0.99,ell=2..10 only',
        'Does not validate the full metric, spin-1/spin-2 sector split, or environmental flux'])
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/paper_kappa_table_audit.json'
    out.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
