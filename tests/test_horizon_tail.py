import sys
from pathlib import Path
import numpy as np
from scipy.integrate import quad_vec
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from environment_horizon_tail import horizon_tail


def test_tail_matches_independent_improper_integral_for_both_frequency_signs():
    for gamma in (.88,-.4):
        cutoff=.003;rp=1.5
        coefficients=np.array([.7+.2j,-.1+.05j,.3-.4j,.03+.08j])
        powers=np.array([0,1,-2j*gamma,1-2j*gamma])
        x=np.geomspace(1.01,9.,30)
        values=np.exp(np.log(x)[:,None]*powers)@coefficients
        correction,report=horizon_tail(rp+cutoff*x,values,rp,gamma,cutoff)
        exact,_=quad_vec(lambda t:cutoff*np.sum(coefficients*np.exp(-(powers+1)*t)),0,np.inf,epsabs=1e-13)
        np.testing.assert_allclose(correction,exact,rtol=2e-10,atol=2e-13)
        assert report['held_out_integrand_residual']<1e-11
