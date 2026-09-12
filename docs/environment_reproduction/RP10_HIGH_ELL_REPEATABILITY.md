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
