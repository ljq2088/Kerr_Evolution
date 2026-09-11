"""Audit continuous radial fields against saved finite Gauss responses."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from environment_response import SampledResponse
from environment_cloud import mode_flux


def pair(z):return [float(z.real),float(z.imag)]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs',nargs='+',type=Path)
    args=parser.parse_args()
    for path in args.inputs:
        data=json.loads(path.read_text())
        response=SampledResponse.from_report(data)
        p=data['parameters'];metric=p['metric']
        g=response.green
        r=np.unique(np.r_[np.linspace(max(g.rmin,2.),150.,601),response.panels])
        field,derivative=response.evaluate(r)
        boundaries=data['radial_response_at_panel_boundaries']
        probe=np.array([row['r'] for row in boundaries])
        old=np.array([complex(*row['field']) for row in boundaries])
        new=response.evaluate(probe)[0]
        omega_c=p['omega']-metric['m_g']/(metric['orbital_radius']**1.5+metric['a'])
        flux=mode_flux(p['omega'],p['scalar_m'],omega_c,1,p['alpha'],metric['a'],
            response.up_coefficient,response.horizon_coefficient)
        report=dict(status='continuous_finite_source_field_not_converged',
            input=path.name,input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            parameters=p,method='Panelwise Legendre source interpolation and bidirectional adaptive Green quadrature',
            up_coefficient=pair(response.up_coefficient),z_h=pair(response.horizon_coefficient),flux=flux,
            original_gauss_flux=data['flux'],
            boundary_field_relative_max_error=float(np.max(abs(new-old))/max(np.max(abs(old)),1e-300)),
            radial_field=[dict(r=float(x),field=pair(f),derivative=pair(d)) for x,f,d in zip(r,field,derivative)],
            limitations=['Finite source cutoffs unchanged; no horizon extrapolation applied',
                'Interpolation needs independent source-grid refinement',
                'One scalar spheroidal multipole, not the full wake'])
        output=path.with_name('continuous_'+path.name)
        output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({k:report[k] for k in ('input','boundary_field_relative_max_error','flux')}),flush=True)


if __name__=='__main__':main()
