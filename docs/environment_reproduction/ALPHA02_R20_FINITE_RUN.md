# alpha=0.2 Kerr r0=20M finite flux run

The full article objective remains all figures, backgrounds, scans and convergence. This queue starts the first L18/nt18 production-resolution alpha=0.2 Kerr orbit; it is not the full Figure 2 or Figure 3 scan.

Requested finite ranges are scalar ell<=6 for propagating infinity flux and ell<=5 for nonzero orbital-effective horizon flux. Their union has 21 independent (ell,m) channels: nine infinity and eighteen horizon channels, with overlap. The m=1 effective flux vanishes and is excluded from these sums. Metric sectors +1 through +5 include their distinct conjugate-driven scalar channels; metric -6 supplies scalar (5,-5) without an unnecessary scalar (7,7) channel.

The preflight alpha02_r20_finite_boundary_preflight.json verifies the six-versus-five-term inverse-r boundary guard at 1000M for every required scalar channel, maximum ratio 1.67e-7. This is not a rigorous boundary error bound.

Parameters: alpha=0.2, threshold Kerr |211> cloud normalized to unit Killing mass, r0=20M, metric reconstruction index through 18, nt18, nr8, h32 logarithmic first panel, inner source offset 0.0005M, source cutoff 320M, Green offset 0.0001M, outer Green radius 1000M, series boundary. Two workers precompute resumable metric samples. Source-domain and near-horizon convergence must be checked for this alpha; results from alpha=0.3 do not establish them.

The sequential queue logs are outputs/alpha02_r20_finite_mg1.log through mg5.log and mg-6.log. Any process error stops later sectors. At launch the corresponding coverage inventory has zero completed production channels and null totals. Completion must be established from terminal batch reports and the background-specific coverage tool, not from this queue declaration.


The first metric-m1 batch scalar_batch_mg1_a095a9ac6eae.json is completed:
scalar m2/ell2,4,6 and m0/ell0,2,4, each with 88 source nodes. Finite source
values, scalar identifiers and equality of response/batch fluxes were checked.
Fingerprints and fluxes are in alpha02_r20_mg1_completed_audit.json.
The finite inventory now has 3/9 infinity and 5/18 horizon channels; totals
remain null. The queue has advanced to metric m2. This is not yet a full
orbit comparison with the paper and does not validate alpha0.2 convergence.
