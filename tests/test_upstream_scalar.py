import sys
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from kerr_scalar import Config, KerrScalar, sin_squared_matrix
from scalar_upstream import convert
import json

# Independent transcription of RIPLEY_OPERATOR.md, set field spin s=0.
def test_coordinate_export_roundtrip(tmp_path):
    mass, chi, length = 2.3, .7, 3.1
    rp = mass*(1+np.sqrt(1-chi**2))
    sigma = np.array([0., .2, 1.])
    tau = np.array([0., 1.])
    u = np.arange(6).reshape(2, 3, 1)*(1+2j)
    source, target = tmp_path/'old.npz', tmp_path/'new.npz'
    np.savez(source, config=json.dumps(dict(mass=mass, spin=chi)),
             sigma=sigma, tau=tau, u=u, p=2*u, ell=[2])
    convert(source, target, length)
    with np.load(target, allow_pickle=False) as result:
        np.testing.assert_allclose(result['T']+4*mass*np.log(rp/mass), tau)
        np.testing.assert_allclose(result['T_elapsed'], tau)
        np.testing.assert_allclose(result['R'], length**2/rp*sigma)
        np.testing.assert_allclose(result['psi0']*length**2, u)
        np.testing.assert_allclose(result['psi0_T']*length**2, 2*u)
        # Away from scri, both conventions reconstruct the same physical field.
        np.testing.assert_allclose(result['R'][None,1:,None]*result['psi0'][:,1:],
                                   sigma[None,1:,None]/rp*u[:,1:])


def test_scalar_matches_upstream_general_spin_limit():
    maximum = 0.
    for M in (1., 2.3):
        for chi in (0., .7, -.8, .99):
            for m in (0, 1, 2):
                obj = KerrScalar(Config(mass=M, spin=chi, m=m, lmax=5, n=24))
                a, rp, Lg = M*chi, obj.horizon, 1.4*M
                k = Lg**2/rp
                R = k*obj.s
                S = sin_squared_matrix(m,5)
                I = np.eye(len(obj.ls))
                ar = 8*M*(2*M-a*a*R/Lg**2)*(1+2*M*R/Lg**2)
                A = ar[:,None,None]*I-a*a*S
                C = -2*(Lg**2-(8*M*M-a*a)*R**2/Lg**2+4*a*a*M*R**3/Lg**4)/k
                D = -(Lg**2-2*M*R+a*a*R**2/Lg**2)*R**2/Lg**2/k**2
                E = 2*(2*M*(4*M*R/Lg**2-3*a*a*R**2/Lg**4)-a*a*R/Lg**2)+2j*a*m*(1+4*M*R/Lg**2)
                F = (2*R*(-1+3*M*R/Lg**2-2*a*a*R**2/Lg**4)+2j*a*m*R**2/Lg**2)/k
                V = (2*(M*R/Lg**2-a*a*R**2/Lg**4)+2j*a*m*R/Lg**2)[:,None]+obj.ls*(obj.ls+1)
                for x,y in [(A,np.linalg.inv(obj.ainv)),(C,obj.C[:,0]),(D,obj.D[:,0]),(E,obj.E[:,0]),(F,obj.F[:,0]),(V,obj.V)]:
                    error = np.max(np.abs(x+y)/(1+np.abs(y)))
                    maximum=max(maximum,float(error))
                    assert error<1e-12
