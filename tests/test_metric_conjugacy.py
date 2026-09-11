import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from environment_lorenz_mode import LorenzMetricMode,ConjugateMetricMode


def test_conjugate_adapter_matches_independent_opposite_frequency_metric():
    for a,r0,r in ((.6,6.,4.5),(.8771530275949366,20.,3.)):
        negative=LorenzMetricMode(r0,a,-1,2)
        positive=LorenzMetricMode(r0,a,1,2)
        reused=ConjugateMetricMode(negative)
        direct=positive(r,1.1)
        error=np.max(abs(reused(r,1.1)-direct))
        assert error<2e-9*np.max(abs(direct))+1e-12
        assert reused.provenance==positive.provenance
        assert reused.omega==positive.omega
