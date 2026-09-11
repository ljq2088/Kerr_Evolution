import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from environment_bound_spectrum import complex_angular_eigenvalue,quasibound_mode
from environment_schwarzschild_cloud import quasibound_211


def test_complex_angular_branch_and_small_parameter_limit():
    c2=1e-7+2e-7j
    lam=complex_angular_eigenvalue(0,0,c2)
    assert abs(lam+c2/3)<1e-14
    np.testing.assert_allclose(complex_angular_eigenvalue(0,0,c2.conjugate()),lam.conjugate(),atol=1e-14)
    assert abs(complex_angular_eigenvalue(2,2,0.)-6)<1e-14


def test_complex_kerr_shooting_recovers_schwarzschild_special_case():
    reference=complex(*quasibound_211(.3)['omega'])
    actual=complex(*quasibound_mode(0.,.3,1,1,seed=reference,rtol=2e-12,outer_efolds=45.)['omega'])
    assert abs(actual-reference)<2e-11


def test_threshold_cloud_is_a_stationary_complex_frequency_root():
    reference=.29629324847975713
    result=quasibound_mode(.8771530275949366,.3,1,1,seed=reference+0j,
                          outer_efolds=45.,offset=1e-6,rtol=2e-12)
    actual=complex(*result['omega'])
    assert abs(actual.real-reference)<2e-11
    assert abs(actual.imag)<2e-12
