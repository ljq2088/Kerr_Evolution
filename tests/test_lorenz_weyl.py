import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lorenz_weyl import integrated_weyl_amplitudes


def test_two_weyl_fields_against_public_table():
    # 2406.12510v3 amplitudes.dat; explicit time-integral convention.
    result=integrated_weyl_amplitudes(4.)
    np.testing.assert_allclose(result[-2],[-.033785123703980326-.14493836520836942j,
                                         .0549449433759872-.13166240267951732j],rtol=2e-10)
    np.testing.assert_allclose(result[2],[-91.99096507683075-238.27417235428504j,
                                        -.18534505760677758-.0962307192687292j],rtol=2e-10)


def test_static_inverse_time_derivative_rejected():
    with pytest.raises(ValueError,match='Static'):
        integrated_weyl_amplitudes(6.,m=0)
