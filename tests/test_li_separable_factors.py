import mpmath as mp
import pytest
from li_separable_factors import sigma_fourier_coefficients,gamma_product_terms,evaluate_terms

@pytest.mark.parametrize('beta',[0,1,2,3])
def test_coefficients_against_independent_integrals(beta):
    with mp.workdps(70):
        r=mp.mpf('1.48');a=mp.mpf('.88')
        f=sigma_fourier_coefficients(str(r),str(a),beta)
        for p in (0,2,6,12):
            expected=(1 if p==0 else 2)/mp.pi*mp.quad(
                lambda t:(r*r+a*a*mp.cos(t)**2)**(-beta)*mp.cos(p*t),[0,mp.pi/2,mp.pi])
            assert abs(f[p]-expected)<mp.mpf('1e-58')
        assert all(f[p]==0 for p in range(1,13,2))

@pytest.mark.parametrize('beta,sigma',[(3,0),(0,3),(3,1),(1,3),(2,2),(0,0)])
def test_gamma_product_and_binomial_factor(beta,sigma):
    with mp.workdps(70):
        r=mp.mpf('1.48');a=mp.mpf('.88');t=mp.mpf('1.13')
        exact=(r+1j*a*mp.cos(t))**(-beta)*(r-1j*a*mp.cos(t))**(-sigma)
        coarse=evaluate_terms(gamma_product_terms(str(r),str(a),beta,sigma),str(t))
        fine=evaluate_terms(gamma_product_terms(str(r),str(a),beta,sigma,pmax=32),str(t))
        assert abs(coarse/exact-1)<mp.mpf('3e-6')
        assert abs(fine/exact-1)<mp.mpf('1e-15')

@pytest.mark.parametrize('beta',[-1,4,1.5,True])
def test_invalid_exponents_rejected(beta):
    with pytest.raises(ValueError):sigma_fourier_coefficients(2,.88,beta)
