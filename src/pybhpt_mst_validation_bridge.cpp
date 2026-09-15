// Isolated diagnostic bridge; never loaded by the production Python process.
#include "mst.hpp"
#include <cmath>
#include <exception>
extern "C" int patched_mst_values(double a, int s, int ell, int m, double omega,
    double lambda, int nr, const double *r, double *meta, double *values) {
    try {
        MstParameters params(a, s, ell, m, 2.0 * omega, lambda);
        nu_solver(params);
        Complex nu = params.getRenormalizedAngularMomentum();
        meta[0] = nu.real(); meta[1] = nu.imag();
        Complex beta = betaMST(0, params), up = alphaRn_cf(1, params), down = gammaLn_cf(-1, params);
        Complex residual = beta + up + down;
        meta[2] = std::abs(residual);
        meta[3] = std::abs(residual)/(std::abs(beta)+std::abs(up)+std::abs(down));
        if(!std::isfinite(nu.real()) || !std::isfinite(nu.imag()) || std::abs(nu) == 0.) return 2;
        MstSeriesWorkspace workspace(params);
        for(int k=0; k<nr; ++k) {
            for(int b=0;b<2;++b) {
                BoundaryCondition bc = b == 0 ? In : Up;
                Complex v = workspace.getNormalizedSolution(bc,r[k]).getValue();
                Complex d = workspace.getDerivativeOfNormalizedSolution(bc,r[k]).getValue();
                values[k*8+b*4] = v.real(); values[k*8+b*4+1] = v.imag();
                values[k*8+b*4+2] = d.real(); values[k*8+b*4+3] = d.imag();
            }
        }
        return 0;
    } catch(const std::exception&) { return 3; }
}
