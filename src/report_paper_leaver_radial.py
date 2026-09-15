"""Cross-check actual scalar In solutions and horizon normalization by series."""
import json,hashlib
from pathlib import Path
import numpy as np
from paper_leaver_radial import LeaverIngoing
from environment_radial import RadialGreen


def main():
    root=Path(__file__).resolve().parents[1];folder=root/'docs/environment_reproduction';rows=[]
    names=['forced_mode_nr8_nt18_L18_mg-1_sl0_inner0.0005_outer320_log_h32.json',
           'forced_mode_nr8_nt18_L18_mg1_sl2_inner0.0005_outer320_log_h32.json',
           'forced_mode_nr8_nt18_L18_alpha0.2_rp20_mg-1_sl0_inner0.0005_outer320_log_h32.json',
           'forced_mode_nr8_nt18_L18_schwarzschild_frozen_mg-1_sl0_inner0.0005_outer320_log_h32.json']
    for name in names:
        data=json.loads((folder/name).read_text());p=data['parameters'];a=p['metric']['a'];mu=p['alpha'];w=p['omega'];ell=p['scalar_ell'];m=p['scalar_m']
        old=RadialGreen(a,mu,w,ell,m,rmax=1000.,offset=1e-4,rtol=1e-11)
        cases=[];expected_current=4*old.rp*(w-m*a/(2*old.rp))
        for terms in (200,800,3200):
            new=LeaverIngoing(a,mu,w,ell,m,terms=terms);samples=[]
            for r in (old.rp+.001,3.,10.,20.,40.):
                actual=new.state(r);reference=old.insol.sol(r)
                current=-2*(r-old.rp)*(r-(2-old.rp))*np.imag(actual[0].conjugate()*actual[1])
                samples.append(dict(r=r,relative_state_max=float(np.max(abs(actual/reference-1))),
                  flux_current_relative=float(current/expected_current-1)))
            cases.append(dict(terms=terms,samples=samples));print(name,terms,max(s['relative_state_max'] for s in samples),flush=True)
        rows.append(dict(input=name,input_sha256=hashlib.sha256((folder/name).read_bytes()).hexdigest(),
          alpha=mu,a=a,omega=w,ell=ell,m=m,expected_inward_charge_current=expected_current,cases=cases))
    result=dict(status='independent_In_series_and_horizon_current_audit',rows=rows,
      limitations=['Checks In solution, not independently the Up solution or full forced Green response',
        'Series convergence slows near x=1; low term counts are retained as unconverged evidence',
        'Analytic unit horizon normalization; no fit to radial ODE data',
        'Wronskian current evaluated in double after high-precision series summation'])
    (folder/'leaver_radial_full_audit.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
