import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from environment_radial import RadialGreen


def test_coulomb_matches_distant_wave_series_normalization():
    # At this nonthreshold frequency the inverse-r series is independently
    # well resolved. Compare normalized solutions after propagation inward.
    reference=RadialGreen(.6,.3,.4,2,2,rmax=8000.,rtol=2e-11)
    coarse=RadialGreen(.6,.3,.4,2,2,rmax=1000.,rtol=2e-11,infinity_method='coulomb')
    fine=RadialGreen(.6,.3,.4,2,2,rmax=4000.,rtol=2e-11,infinity_method='coulomb')
    radii=np.array([8.,20.,50.])
    target=reference.upsol.sol(radii)[0]
    error=lambda g:np.max(abs(g.upsol.sol(radii)[0]/target-1))
    assert error(fine)<error(coarse)/8
    assert error(fine)<2e-6


def test_coulomb_signed_frequency_reality():
    positive=RadialGreen(.6,.3,.4,2,2,rmax=1000.,infinity_method='coulomb')
    negative=RadialGreen(.6,.3,-.4,2,-2,rmax=1000.,infinity_method='coulomb')
    r=np.array([5.,20.,100.])
    np.testing.assert_allclose(negative.upsol.sol(r),positive.upsol.sol(r).conjugate(),rtol=1e-10,atol=1e-12)
