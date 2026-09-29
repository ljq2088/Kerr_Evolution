# Fig.10 follow-up, 2026-09-28

## Execution
The Fig.10 orbital scan has been prioritized over remaining Fig.9 high-m modes. Completed outputs and metric cache samples were retained. Current process status is execution.json. The first new radius is 18.3M, where the paper reports a horizon-mode dominance exchange.

## Stronger normalization hypothesis
Define the unshifted Leaver radial series by
R_raw = exp(q*r)*(r-r_minus)^beta*x^s*sum(a_n*x^n), a0=1,
x=(r-r_plus)/(r-r_minus).

The local normalized cloud is
R_norm = A*exp(q*(r-r_plus))*((r-r_minus)/(r_plus-r_minus))^beta*x^s*sum(a_n*x^n),
where A=1/sqrt(M_raw_shifted) normalizes the finite BL cloud energy.

Compared with the physical flux in epsilon^2*zeta^2 units, the independently predicted raw-series flux multiplier is
K = alpha^6*M_raw_shifted*abs(exp(q*r_plus)*(r_plus-r_minus)^beta)^2.

This yields K_211=463.91056327644213 and K_322=105414901.77833593. They are computed from the background solution, not fitted to Fig.10. At r0=20M, all eight plotted infinity channels then agree with PDF readouts within 0.20%--1.97%; the cloud211 subset has 0.30%--0.81% residuals. This explains both enormous scale offsets much better than a universal (r0*u^t)^2. K is independent of orbital radius, so the new 18.3M calculation is an important falsifiable check.

This remains a hypothesis about the published plotting convention, not confirmation of unpublished author code. Physical source, response, and flux outputs remain unchanged. Horizon modes do NOT all agree after this transformation.

## Independent controls completed
1. Both cloud radial profiles match Fig.3 in the declared mass normalization (~1e-5 relative vector-readout error).
2. Selected cloud322 separated sources agree with direct covariant contraction to 1e-14 for the same metric inputs.
3. Fixed-source Green boundary changes: outer 500,1000,1500M and horizon offset 1e-4 vs1e-5M change the tested cloud322 horizon fluxes by less than 7e-11 relatively. This rules out these Green endpoint choices as the cause of the large residual. It does not test source resolution/metric accuracy.
4. Using only the printed Appendix terms instead of the covariant-identity-corrected source does not restore the missing horizon amplitudes; it is retained as a diagnostic, never adopted as a correction.
5. The available cloud322 infinity modes at r20 sum to 0.0005664214845, while Fig.13 gives approximately 0.0005663725626 (0.00864% difference). This is a partial sum of dominant modes and a digitized total, not a completed full-mode total comparison. Cloud1 Fig.13 curves are clipped at this radius and cannot be numerically read from the clipped path.

## Open problems
Horizon-channel discrepancies, metric/source quadrature convergence, source near-horizon roundoff and full orbit dependence remain unresolved. Do not claim full Fig.10 reproduction. The raw-Leaver comparison is figure10_raw_leaver_control.png, with exact factors and pointwise comparisons in raw_leaver_normalization_control.json.
