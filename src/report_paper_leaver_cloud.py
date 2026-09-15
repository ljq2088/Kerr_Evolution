"""Independent cloud spectrum, radial profile and normalization audit."""
import json,hashlib
from pathlib import Path
import mpmath as mp
import numpy as np
from paper_leaver_cloud import LeaverThresholdCloud
from environment_source import ThresholdCloud


def main():
    rows=[]
    for alpha in (.2,.3):
        production=ThresholdCloud(alpha=alpha);cases=[]
        for terms in (100,200,400):
            new=LeaverThresholdCloud(alpha,terms=terms)
            e,q=new.integrals(192);coarse,_=new.integrals(96)
            radii=[production.rp+.001,3.,10.,20.,40.,100.,200.,320.]
            samples=[]
            for radius in radii:
                old=production.radial(radius)[:,0];current=new.radial(radius)/np.sqrt(e)
                samples.append(dict(r=radius,leaver=current.tolist(),shooting=old.tolist(),
                    field_relative=float(current[0]/old[0]-1),derivative_relative=float(current[1]/old[1]-1)))
            cases.append(dict(terms=terms,a=mp.nstr(new.a,40),omega=mp.nstr(new.omega,40),
                continued_fraction_residual=mp.nstr(new.residual,8),
                energy_charge_relative=float(e/(float(new.omega)*q)-1),
                energy_quadrature_relative=float(coarse/e-1),samples=samples))
            print(alpha,terms,'a error',float(new.a)-production.a,'R max',max(abs(s['field_relative']) for s in samples),flush=True)
        rows.append(dict(alpha=alpha,production_a=production.a,production_omega=production.omega,cases=cases))
    root=Path(__file__).resolve().parents[1]
    result=dict(status='independent_threshold_Leaver_cloud_comparison',rows=rows,
      implementation_sha256={f:hashlib.sha256((root/'src'/f).read_bytes()).hexdigest() for f in ['paper_leaver_cloud.py','environment_cloud.py','environment_source.py']},
      limitations=['Real synchronous |211> mode only',
        'Angular eigensystem and recurrence use 50 decimal digits; stress-energy quadrature is double',
        'Comparison uses independent unit Killing-mass normalization, no fitted amplitude',
        'This checks the cloud input, not the forced metric/scalar response'])
    (root/'docs/environment_reproduction/leaver_cloud_full_audit.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
