import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from environment_lorenz_mode import LorenzMetricMode,ConjugateMetricMode
from environment_metric_sampling import precompute_metric


def test_parallel_tensor_samples_match_direct_and_resume(tmp_path):
    base=LorenzMetricMode(6.,.6,2,2)
    radii=[4.5,9.];theta=[.7,1.2]
    sampled,audit=precompute_metric(base,radii,theta,tmp_path,workers=2)
    assert audit['reconstructed_radii']==2
    for r in radii:
        for t in theta:
            expected=base(r,t)
            # Independent backend/auxiliary-solver cache histories can differ
            # at ~5e-11 relative; compare within reconstruction accuracy, not
            # bitwise arithmetic. Resume below still requires exact equality.
            np.testing.assert_allclose(sampled(r,t),expected,rtol=1e-9,atol=1e-12)
            np.testing.assert_allclose(ConjugateMetricMode(sampled)(r,t),expected.conj(),rtol=1e-9,atol=1e-12)
    again,resumed=precompute_metric(base,radii,theta,tmp_path,workers=2)
    assert resumed['cached_radii']==2 and resumed['reconstructed_radii']==0
    assert again.provenance==base.provenance
    # Modifying a returned tensor must not mutate the durable or memory cache.
    again(4.5,.7)[:]=0
    np.testing.assert_array_equal(again(4.5,.7),sampled(4.5,.7))
    file=next(tmp_path.rglob('*.npz'))
    with file.open('wb') as stream:
        np.savez(stream,metadata='wrong',r=4.5,h=np.zeros((2,4,4)))
    with pytest.raises(ValueError,match='provenance'):
        precompute_metric(base,radii,theta,tmp_path,workers=2)
