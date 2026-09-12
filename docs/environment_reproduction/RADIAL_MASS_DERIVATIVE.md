# Auxiliary mass-squared radial variation (diagnostic)

At fixed real frequency, write nu for the auxiliary mass squared and A(nu) for
its scalar angular eigenvalue. The homogeneous radial equation is

    Delta R'' + Delta' R' + V R = 0,
    V = K^2/Delta - nu r^2 - a^2 omega^2 + 2 a m omega - A(nu).

Differentiation with U = partial_nu R gives

    Delta U'' + Delta' U' + V U = (r^2 + A_nu) R.

`environment_radial_variation.py` integrates both equations together. Forward
derivative arithmetic differentiates all coefficients of the four-correction
horizon Frobenius expansion and the six-correction outgoing expansion. In the
latter, k_nu=-1/(2k) and p=-1+i(2 omega^2-nu)/k; normalization phases are also
differentiated. The horizon exponent is independent of nu at fixed omega.
The angular derivative comes from the independently tested eigenpair variation.

For W=Delta(R_in R_up' - R_up R_in'), the differentiated product gives W_nu.
Both W and W_nu must be independent of radius. Three checks use a=0.6 or 0.877,
omega=0.1, 0.03 or -0.1, ell=1,3,2, and m=1,1,-1. On 40 radii to 300M,
comparison with the independent production solver and fourth-order centered
mass differences gave original-solution relative discrepancies below 4e-12,
tangent discrepancies below 4.4e-9, and W/W_nu spatial spreads below 1.9e-11.
These are numerical checks at the tested finite boundaries, not outer-boundary
convergence estimates or a proof for all multipoles.

Scope: real nonzero propagating auxiliary frequencies, subextremal Kerr.
No bound/Coulomb derivative is implemented. This module is not imported into
production Lorenz reconstruction. Differentiating the forced trace source,
Green coefficients and full separated Taylor jets is still required before
this can replace the finite-difference kappa construction. These checks do not
resolve the disagreement with the paper's fluxes.
