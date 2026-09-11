import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from environment_schwarzschild_cloud import SchwarzschildCloud


def test_complex_cloud_kg_and_horizon_decay_balance():
    for alpha in (.2,.3):
        cloud=SchwarzschildCloud(alpha=alpha)
        balance=cloud.decay_balance()
        assert cloud.spectral_omega.imag<0
        assert abs(balance['energy']-1)<1e-10
        assert abs(balance['energy_balance']/balance['horizon_energy'])<2e-6
        assert abs(balance['charge_balance']/balance['horizon_charge'])<2e-6
        for r in (3.,10.,60.):
            field,H,inverse=cloud.hessian(r,1.1)
            residual=np.einsum('ij,ij->',inverse,H)-alpha**2*field
            assert abs(residual/field)<1e-11


def test_freezing_time_frequency_has_explicit_kg_defect():
    cloud=SchwarzschildCloud(alpha=.3,freeze_decay=True)
    for r in (3.,10.,60.):
        field,H,inverse=cloud.hessian(r,1.1)
        residual=np.einsum('ij,ij->',inverse,H)-cloud.mu**2*field
        expected=(cloud.omega**2-cloud.spectral_omega**2)/(1-2/r)*field
        np.testing.assert_allclose(residual,expected,rtol=2e-8,atol=1e-13)
        assert abs(expected)>1e-10
