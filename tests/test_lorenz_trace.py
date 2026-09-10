import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lorenz_trace import trace_amplitudes


def test_trace_against_published_complex_amplitudes():
    # 2406.12510v3 source archive, amplitudes.dat, r0=6, a=.6, ell=m=2.
    reference=np.array([-.4980022825100904+.2753958041341394j,
                         .010397600098005078+.013217659834909017j])
    baseline=np.array(trace_amplitudes(6.,rmax=1000.,offset=1e-5))
    refined=np.array(trace_amplitudes(6.,rmax=2000.,offset=1e-6))
    np.testing.assert_allclose(refined,reference,rtol=2e-8,atol=1e-10)
    assert np.max(abs(refined-baseline)) < 1e-8
