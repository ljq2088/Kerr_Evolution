# High-L spin-2 repeatability diagnostic

The sector projection experiment re-evaluated L17 and L18 contributions at the
same q18/r0=10/a=.8771530275949366/m1 settings, using analytic kappa. Its sector
sum does not reproduce the corresponding increment of the earlier full run:
maximum derivative differences are about 1.24e-4 and 1.18e-4. The sector data
must therefore not be treated as a reliable decomposition of that prior run.
Within this re-evaluation the spin-2 sector has much larger projected jumps
than spin-1/0, but that observation alone does not locate the physical error.

Three calls to spin2_metric with identical inputs r=10.00005, theta=1.1,
order8 exhibit nonzero differences at L17 and L18 under the default angular
backend. The analogous experiment in a fresh process using the existing dense
real angular diagnostic gives exactly identical arrays on all three calls.
The radial data and source amplitudes are cached inside each process. This
provides evidence that angular evaluation repeatability contributes to the
observed high-L variability. It is not proof that angular error explains the
paper flux discrepancy, nor does exact repeatability prove angular accuracy.
The default and dense repeat measurements are saved in the two accompanying
rp10_*_spin2_repeatability.json files.

A full L1..18 q18 matching comparison with dense angular functions and analytic
kappa is the next required diagnostic. Production reconstruction is unchanged.


## Accuracy of the dense angular functions at the tested high-L parameters

At the actual r0=10M, a=.8771530275949366, m=1 frequency, L17 and L18
were tested for spin weights -2,-1,0,1,2 on 100 Gauss-Legendre nodes.
Each dense eigenproblem was also enlarged by 20 basis functions, using
the same symmetric eigensolver and phase convention. The maximum relative
change in S, its first derivative and its second derivative was
1.465514541e-15. The normalized angular ODE residual was at most
2.537246019e-14, and unit-sphere normalization error at most 7.793765633e-14.
The separation eigenvalue differed from the independent radial C++ interface
by at most 4.547473509e-13. All ten cases and dimensions are retained in
rp10_high_ell_dense_angular_accuracy.json.

These tests support angular basis convergence for these functions at the tested
parameters. They do not bound cancellation amplification through metric
reconstruction, radial errors or environmental fluxes. The full dense-angular
matching sequence remains a separate running calculation.


## Complete dense-angular matching sequence

The q18 L1..18 calculation with analytic kappa and dense angular functions
has completed. Physical, projection and Taylor parameters match the preceding
analytic-kappa run; only the angular backend changes. At L18, the value jump
is 6.3638530857e-06, compared with 0.000153986363168;
the derivative mismatch is 1.76897217777e-05, compared with
0.000108257434104. Ratios dense/default are
0.0413273809 and 0.163404222.

The paired sequence and input hashes are in rp10_metric_matching_dense_angular_audit.json.
These are projected metric residuals; no inference about scalar-flux changes
is justified without recomputing the source. Nonzero high-order residuals still
require quadrature, arithmetic and radial-input checks. Production remains unchanged.


## Dense angular functions with extended Taylor arithmetic

The complete L1..18 q18 run is now available with both dense angular functions
and extended Taylor coefficients, keeping analytic kappa. At L18, the value
jump is 4.37068903337e-06 (double 6.3638530857e-06),
and derivative mismatch is 1.6890593223e-05
(double 1.76897217777e-05). Their extended/double ratios are
0.686799173 and 0.954825262.
Radial and angular inputs remain double in both runs.

Input fingerprints and all 18 paired rows: rp10_metric_matching_dense_precision_audit.json.
This comparison does not establish quadrature convergence, radial accuracy,
or scalar flux changes. The independent dense-source flux run is pending.
