import sys
from pathlib import Path
import numpy as np
from scipy.signal import convolve2d
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lorenz_jet import Jet


def test_triangular_product_matches_independent_full_convolution():
    rng=np.random.default_rng(704)
    for order in (0,4,6,8,10):
        shape=(order+1,order+1)
        a=Jet(order=order,coefficients=rng.normal(size=shape)+1j*rng.normal(size=shape))
        b=Jet(order=order,coefficients=rng.normal(size=shape)+1j*rng.normal(size=shape))
        reference=convolve2d(a.c,b.c)[:order+1,:order+1]
        for i in range(order+1):
            for j in range(order+1):
                if i+j>order:
                    reference[i,j]=0
        np.testing.assert_allclose((a*b).c,reference,rtol=1e-13,atol=1e-13)
