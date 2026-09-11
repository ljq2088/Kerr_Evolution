"""Isolate homogeneous angular re-solve noise in high-L Lorenz tensors.

Diagnostic only: production reconstruction and running metric jobs are unchanged.
Source amplitudes already cached by the baseline are deliberately held fixed.
"""
import json
from pathlib import Path
import numpy as np
from pybhpt.radial import RadialTeukolsky
from pybhpt.swsh import SpinWeightedSpheroidalHarmonic
from environment_angular_diagnostic import DenseRealHarmonic
import lorenz_metric
import lorenz_chi
import lorenz_spin1
import lorenz_spin1_chiral


def main():
    a=.8771530275949366;r0=20.;r=30.;theta=1.1;m=2
    omega=m/(r0**1.5+a)
    angular_checks=[]
    nodes,weights=np.polynomial.legendre.leggauss(80)
    angles=np.arccos(nodes)
    for ell in (2,6,18):
        for spin in (-2,-1,0,1,2):
            harmonic=DenseRealHarmonic(spin,ell,m,a*omega)
            radial=RadialTeukolsky(spin,ell,m,a,omega,np.array([r]))
            value=harmonic(angles);derivative=harmonic(angles,deriv=1)
            second=harmonic(angles,deriv=2)
            # pybhpt reports lambda=A+gamma^2-2*m*gamma.
            gamma=a*omega;A=harmonic.eigenvalue-gamma**2+2*m*gamma
            potential=gamma**2*nodes**2-2*gamma*spin*nodes+spin+A-(m+spin*nodes)**2/(1-nodes**2)
            residual=second+nodes/np.sqrt(1-nodes**2)*derivative+potential*value
            angular_checks.append(dict(ell=ell,spin=spin,
                cpp_radial_eigenvalue_difference=float(harmonic.eigenvalue-radial.eigenvalue),
                norm=float(2*np.pi*np.dot(weights,abs(value)**2)),
                angular_ode_scaled_residual=float(np.max(abs(residual))/max(1.,np.max(abs(potential*value))))))
    modules=(lorenz_metric,lorenz_chi,lorenz_spin1,lorenz_spin1_chiral)
    rows=[]
    for ell in (2,6,18):
        def field():
            _,h=lorenz_metric.nonstatic_metric(r,theta,r0,a,ell,m,order=6)
            return np.array([[v.value for v in row] for row in h])
        old=np.array([field() for _ in range(4)])
        try:
            for module in modules:module.SpinWeightedSpheroidalHarmonic=DenseRealHarmonic
            new=np.array([field() for _ in range(4)])
        finally:
            for module in modules:module.SpinWeightedSpheroidalHarmonic=SpinWeightedSpheroidalHarmonic
        scale=float(np.max(abs(old[0])))
        rows.append(dict(ell=ell,
            arpack_repeat_absolute_spread=float(np.max(abs(old-old[0]))),
            arpack_repeat_relative_spread=float(np.max(abs(old-old[0]))/scale),
            dense_repeat_absolute_spread=float(np.max(abs(new-new[0]))),
            dense_relative_shift=float(np.max(abs(new[0]-old[0]))/scale),
            arpack_tensors=np.stack((old.real,old.imag),axis=-1).tolist(),
            dense_tensors=np.stack((new.real,new.imag),axis=-1).tolist()))
        print({k:v for k,v in rows[-1].items() if not k.endswith('tensors')},flush=True)
    result=dict(status='angular_repeatability_diagnostic_not_flux_convergence',
        parameters=dict(a=a,r0=r0,r=r,theta=theta,m=m,taylor_order=6,repeats=4),
        angular_checks=angular_checks,tensor_checks=rows,
        limitations=[
            'Existing baseline source-amplitude caches held fixed; not a fresh full dense solve',
            'Radial inputs remain double precision',
            'Repeatability is not an accuracy bound or a modal flux error estimate',
            'Dense eigensolver is only used inside this diagnostic',
            'ARPACK results may vary between invocations'])
    output=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/angular_repeatability_audit.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(output,flush=True)


if __name__=='__main__':main()
