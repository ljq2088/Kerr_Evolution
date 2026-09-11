"""Independent leading local cusp from the Lorenz point-mass singular field.

The leading spherical-mode prediction is compared diagnostically with the
finite spheroidal-mode data; this is not a proof of their asymptotic equality.
"""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import ellipe,eval_legendre
from environment_source import ThresholdCloud,kerr_metric


def main():
    folder=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    path=folder/'particle_field_L18.json';data=json.loads(path.read_text())
    cloud=ThresholdCloud();r0=data['evaluation']['r'];theta=np.pi/2
    op=1/(r0**1.5+cloud.a);g=kerr_metric(r0,theta,cloud.a)
    ut=1/np.sqrt(-(g[0,0]+2*op*g[0,3]+op*op*g[3,3]))
    u=ut*np.array([1.,0.,0.,op]);uc=g@u
    np.testing.assert_allclose(u@g@u,-1.,rtol=0,atol=1e-13)
    phi,hessian,inverse=cloud.hessian(r0,theta)
    along=np.einsum('a,b,ab->',u,u,hessian)
    direct=-(ut*(cloud.omega-op))**2*phi
    np.testing.assert_allclose(along,direct,rtol=1e-12,atol=1e-15)
    A=cloud.mu**2*phi+2*along
    # Check the algebra of the Lorenz singular-field contraction separately.
    hlocal=2*(g+2*np.outer(uc,uc))
    contraction=np.einsum('ij,ij->',inverse@hlocal@inverse,hessian)
    np.testing.assert_allclose(contraction,2*A,rtol=1e-12,atol=1e-15)
    qa=g[2,2];qb=g[3,3]+uc[3]**2
    qbar=quad(lambda beta:np.sqrt(qa*np.cos(beta)**2+qb*np.sin(beta)**2),0,2*np.pi,
              epsabs=1e-12,epsrel=1e-12)[0]/(2*np.pi)
    major=max(qa,qb);minor=min(qa,qb)
    ellipse=2*np.sqrt(major)*ellipe(1-minor/major)/np.pi
    np.testing.assert_allclose(qbar,ellipse,rtol=1e-13)
    predicted=A*qbar/cloud.mu**3
    rows=[]
    for row in data['multipoles']:
        if not row['complete']:continue
        ell=row['ell'];L=ell+.5
        field=complex(*row['particle_field_sum'])/cloud.mu**3
        coefficient=-L*L*field
        # Exact Legendre coefficient of the illustrative isotropic chord cusp.
        integral=(2*ell+1)/2*quad(lambda x:np.sqrt(2*(1-x))*eval_legendre(ell,x),
            -1,1,epsabs=1e-12,epsrel=1e-12)[0]
        exact=-4/((2*ell-1)*(2*ell+3))
        np.testing.assert_allclose(integral,exact,rtol=1e-9,atol=1e-12)
        rows.append(dict(ell=ell,scaled_coefficient=[coefficient.real,coefficient.imag],
            relative_difference=float(abs(coefficient-predicted)/abs(predicted)),
            chord_model_coefficient=exact))
    result=dict(status='leading_local_cusp_diagnostic_not_full_field_validation',
        input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        parameters=dict(a=cloud.a,alpha=cloud.mu,rp=r0,cloud_mass=1.,omega_c=cloud.omega,
                        omega_p=op,u_t_contravariant=ut),
        background_field=float(phi),u_hessian_u=[along.real,along.imag],
        A=[A.real,A.imag],mean_angular_distance_factor=qbar,
        predicted_spherical_coefficient_per_epsilon_q=[predicted.real,predicted.imag],
        rows=rows,limitations=[
            'Prediction derives from local Lorenz h_ab=2m(g_ab+2u_a u_b)/s',
            'Leading large-ell spherical asymptotic; data are finite spheroidal shell sums',
            'Spheroidicity depends on m; no exact basis-equivalence assertion',
            'Subleading local terms and missing higher reconstructed modes remain',
            'No normalization has been fitted to the paper or numerical shells'])
    (folder/'local_scalar_cusp.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(predicted=result['predicted_spherical_coefficient_per_epsilon_q'],rows=rows)))


if __name__=='__main__':main()
