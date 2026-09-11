"""Reintegrate saved actual sources with independently varied outer boundaries."""
import json
from pathlib import Path
import numpy as np
from environment_radial import RadialGreen
from environment_cloud import mode_flux


def main():
    directory=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    out=directory/'threshold_boundary_refinement.json'
    rows=[]
    for orbit in (41.6,41.8):
        source=directory/f'forced_mode_nr2_nt6_L2_alpha0.3_rp{orbit:g}_mg1_sl2_log_h8.json'
        data=json.loads(source.read_text())
        p=data['parameters']; a=p['metric']['a']; omega=p['omega']
        radii=np.array([s['r'] for s in data['samples']])
        integrand=np.array([s['weight']*complex(*s['source']) for s in data['samples']])
        for method,rmax in [('series',1000.)]+[('coulomb',r) for r in (1000.,2000.,4000.,8000.,16000.,32000.)]:
            g=RadialGreen(a,.3,omega,2,2,rmax=rmax,offset=1e-4,rtol=1e-11,infinity_method=method)
            zi=np.sum(g.insol.sol(radii)[0]*integrand)/g.w0
            zh=np.sum(g.upsol.sol(radii)[0]*integrand)/g.w0
            row=dict(orbit=orbit,method=method,rmax=rmax,
                series_last_term_relative=g.series_last_term_relative,
                z_h=[float(zh.real),float(zh.imag)],
                z_inf=[float(zi.real),float(zi.imag)] if g.propagating else [0.,0.],
                wronskian_spread=float(np.max(abs(g.wronskian(radii)/g.w0-1))),
                flux=mode_flux(omega,2,omega-1/(orbit**1.5+a),1,.3,a,zi,zh))
            rows.append(row)
            out.write_text(json.dumps(dict(status='boundary_diagnostic_low_resolution_source',rows=rows),indent=2)+'\n')
            print(json.dumps(row),flush=True)


if __name__=='__main__':
    main()
