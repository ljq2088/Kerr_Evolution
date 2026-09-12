# r0=10M production-parameter metric source matching

The dominant scalar (0,0) horizon response is driven by metric m_g=-1, related by reality to m_g=+1. The earlier source-junction benchmarks at r0=6 and 41.6 do not establish the matching accuracy at r0=10. An independent projected junction run has been started at a=0.8771530275949366, r0=10, m_g=1, reconstruction index 1..18, with 18 Gauss angular nodes.

The direct point-particle stress tensor fixes the coordinate derivative jump

\[
[\partial_r h_{ab,m}]= -\frac{8}{u^t\Delta_0}
\left(u_a u_b+\frac12 g_{ab}\right)
\frac{\delta(\theta-\pi/2)}{\sin\theta}.
\]

The factor 8 includes the Fourier delta normalization. The target is computed from the background circular four-velocity, independently of the reconstructed metric. Ten symmetric components are projected against P_1^1 and P_2^1, using the globally smooth vector sin(theta) partial_theta for polar tensor slots. Jets of order eight provide values and radial derivatives through second order at r0 +/- 5e-5, extrapolated to common radii. These remain a finite set of smooth weak tests, not a proof for every tensor test function or an error bound on fluxes.

Command: src/report_metric_tensor_matching.py --quadrature 18 --ellmax 18 --m 1 --r0 10 --a .8771530275949366 --smooth.

Log: outputs/metric_matching_r10_m1_L18_q18.log. Partial cases are saved after each reconstructed degree in metric_tensor_matching_q18_smooth_m1_r10_a0.877153027595.json. The run must be confirmed terminal with all eighteen cases before reporting its final outcome. No reconstruction coefficient is adjusted to fit the jump, and no existing flux file is changed.


## Completed double-arithmetic sequence

The process completed all eighteen degrees. At cutoff 18 the maximum extrapolated continuity residual is 1.539863477028e-04, and the derivative-jump error is 1.082572794235e-04. The largest expected derivative-jump component is 5.232365106801e+00, giving a ratio 2.068993222257e-05 to that overall target scale. This is not a per-component relative bound or a scalar-flux error bound.

Both diagnostics improve through cutoff 8 (continuity 2.3423e-7, derivative error 4.0977e-8), then rise. The complete sequence and fingerprint are saved in rp10_metric_matching_L18_audit.json and rp10_metric_matching_L18_sequence.png; the plot was visually checked. The high-cutoff rise prevents a convergence claim.

An otherwise identical --extended-jet run has started. It changes Taylor coefficient arithmetic to numpy.clongdouble, while radial and angular input data remain double precision. This is a diagnostic of arithmetic sensitivity, not a fully high-precision reconstruction. Log: outputs/metric_matching_r10_m1_L18_q18_extended.log. Original flux data and the production reconstruction are unchanged.
