# r0=10M inner source cutoff halving

Completed comparison of six channels at L18, angular order18, radial order8, and horizon order64. The source inner offset changes from 0.0005M to 0.00025M; the Green inner offset remains 0.0001M. All other batch and response parameters agree. The 64 normal-panel radii, weights and complex source values agree exactly in every channel. The logarithmic first panel is remeshed, so cutoff and its quadrature change together.

| (ell,m) | Horizon relative change | Infinity relative change |
|---|---:|---:|
| (0,0) | 0.0001903220884 | zero in both |
| (2,0) | 0.00208328474 | zero in both |
| (2,2) | 0.001110989949 | -2.2056419179429078e-06 |
| (4,0) | 0.008082595623 | zero in both |
| (4,2) | 0.007684753016 | 6.737932203802397e-11 |
| (6,2) | 0.003327823703 | 8.423976809123478e-15 |

The dominant (0,0) horizon flux changes from -6.272275357554517e-5 to -6.273469110099658e-5, an increase in magnitude of 0.0190322%. This finite change is far smaller than the approximately 40.7% current finite-sum/paper difference and increases the discrepancy; it does not establish a bound on all uncomputed near-horizon contributions.

The independent two-phase horizon extrapolation with first-order coefficients and window10 gives -6.275067924549729e-05 at the original cutoff and -6.275061059731141e-05 at the smaller cutoff, relative difference -1.09398315e-06. This is encouraging consistency of extrapolated estimates, not a certified error estimate. At the smaller cutoff, held-out absolute integrand residuals are approximately 1e-4, larger than those in the original-cutoff fit. All fit orders/windows are retained in rp10_scalar00_halved_cutoff_tail_audit.json; the fit does not replace source-resolution validation.

The full six-channel comparison and input fingerprints are in rp10_inner_cutoff_halving.json. No normalization fit or production replacement was made.
