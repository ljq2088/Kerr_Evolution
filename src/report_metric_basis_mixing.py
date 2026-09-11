"""Quantify angular-basis leakage only, not a metric or flux error estimate."""
import json
from pathlib import Path
import numpy as np
from environment_angular_diagnostic import DenseRealHarmonic


def main():
    a=.8771530275949366
    orbit=20.
    rows=[]
    for m in (2,5,7,9,11):
        c=a*m/(orbit**1.5+a)
        for spin in (-2,-1,0,1,2):
            for ell in (12,18):
                mode=DenseRealHarmonic(spin,ell,m,c)
                degrees=mode.lmin+np.arange(len(mode.coeffs))
                norm=float(np.sum(abs(mode.coeffs)**2))
                tail=float(np.sum(abs(mode.coeffs[degrees>18])**2))
                rows.append(dict(spin=spin,ell=ell,m=m,spheroidicity=c,
                    angular_norm=norm,spherical_above_18_norm_squared=tail,
                    spherical_above_18_norm=float(np.sqrt(tail)),
                    largest_off_diagonal_coefficient=float(np.max(abs(mode.coeffs[degrees!=ell])))))
                if abs(norm-1)>1e-12:raise ValueError('Angular coefficient normalization failed')
    result=dict(status='angular_basis_comparison_not_metric_truncation_error',
        a=a,orbital_radius=orbit,rows=rows,
        limitations=['Individual angular factors, before metric reconstruction and source contraction',
            'Missing separated modes above 18 are not computed',
            'Does not determine the authors component basis or explain flux discrepancies'])
    output=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/metric_basis_mixing.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'largest_angular_tail_norm':max(row['spherical_above_18_norm'] for row in rows)}))


if __name__=='__main__':main()
