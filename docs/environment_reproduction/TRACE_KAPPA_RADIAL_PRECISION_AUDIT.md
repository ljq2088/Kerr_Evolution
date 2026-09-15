# Actual s=0 trace/kappa radial audit, 2026-09-15

The actual metric mode is ell=1, m=1, rp=20, a=0.8771530275949366, with omega_g=0.011071760582072043. This is the massless metric trace and its auxiliary mass variation; it is separate from the physical massive environment response at omega_s=omega_cloud+m_g Omega_p.

## Formula and normalization review

`lorenz_kappa.scalar_resolvent` uses the scalar angular eigenvalue A and integrates the radial equation with lambda_radial=A+a^2 omega_g^2-2 a m omega_g. The numerical pybhpt radial eigenvalue agrees with this convention to 1.14e-14. The jump source is -16 pi S(pi/2)/u^t, and In/Up amplitudes use the invariant Green products divided by Delta(In Up'-Up In'). The local Taylor recurrence reconstructs higher derivatives from this same radial lambda; it does not confuse A and lambda.

The construction satisfies (Box-nu) h_nu=16 pi T and kappa=-i omega_g partial_nu h_nu at nu=0. The existing finite difference uses a fourth-order Richardson central derivative with step=0.01 omega_g^2=1.2258388238672426e-6. All auxiliary solves stay on the propagating branch; no mass-threshold crossing occurs.

The independent variational implementation includes derivatives of the angular eigenvalue, angular source projector, output angular function, both radial boundary expansions, amplitudes and Wronskian. The radial variation source is (r^2+partial_nu A)R; omitting these normalization/boundary/projector derivatives would compute a different resolvent derivative. These terms are present. Finite-ell angular-projector variation remains a shared theoretical convention, not independently proved equivalent to every author reconstruction.

## New actual trace comparison

At nine original source radii plus the particle radius, the complex sourced trace and its first derivative are compared to invariant products of independent pybhpt s=0 homogeneous solutions:

| Backend | Maximum relative trace value difference | Maximum relative first derivative difference |
|---|---:|---:|
| AUTO / HBL | 1.54e-11 | 1.81e-11 |
| TEUK | 5.46e-11 | 5.78e-11 |
| Forced MST | nonfinite In and Up output | rejected |

The installed pybhpt forced MST path returned NaN values and derivatives for this s=0 problem, confirmed directly at r=1.480716,3,20,40. It is not a valid comparison or a precision improvement. The production trace path uses its own scalar RadialGreen solver, not this failing forced-MST call. The report detects nonfinite backend values before constructing or accepting a Green function.

## New trace+kappa changes propagated to actual complex source

Nine actual radial nodes are retained, with the original 18-point angular projection. The original complete L1 metric source is recomputed at these nodes; agreement with stored source is within 1.51e-11 relative. Only trace and kappa are then changed, keeping spin1 and compact chi fixed. The opposite metric frequency is obtained by conjugating the metric correction alone, leaving the cloud unchanged.

| Variation | Maximum relative change of projected complex scalar (0,0) source |
|---|---:|
| Halve mass difference step | 2.44e-8 |
| Finite difference: rtol=5e-14, rmax=4000 | 3.05e-8 |
| Analytic variational ODE at original boundaries | 3.08e-8 |
| Analytic variation: rmax=8000, offset=1e-5, horizon order8, infinity order12, rtol=5e-14 | 3.09e-8 |

Local trace and kappa jets at theta=1.1 are also compared through fourth radial derivatives and second mixed/angular derivatives. Maximum relative trace-jet changes are about 1.61e-10 and kappa changes about 8.11e-8. All complex values and derivatives are stored in the JSON. This is a bounded node audit, not a new whole-domain flux integral.

The earlier complete 88x18 kappa-only replacement already yielded a horizon-flux change of about -2.245e-9 across analytic rmax=2000,4000,8000 (`dipole_kappa_flux_boundary_audit.json`). That old result is cited, not repeated or presented as new.

## A real near-horizon precision limit

A separate refinement probe using rmax=8000, endpoint orders4/6, offset=1e-6 and rtol=5e-14 caused the backward variational Up integration to fail at r-r_plus=1.06163e-6 after approximately 1.48 million RHS evaluations: the required step was smaller than floating-point spacing. The In integration completed. Offset=1e-5 at the same tolerance, and offset=1e-4, completed normally.

This demonstrates a limitation of integrating the singular BL radial state too close to the horizon with excessive tolerance demands. It is not the cause of the historical discrepancy, because production uses offset=1e-4 and completes successfully. Further approach to the horizon should factor the Frobenius phase or use a logarithmic/tortoise integration coordinate, instead of reducing the BL endpoint blindly.

No new evidence supports a percent-level s=0 trace/kappa radial integration error at this actual dipole point. Shared physical source normalization, angular truncation and global reconstruction conventions remain outside these numerical comparisons.

Reproduce the completed audit with `PYTHONPATH=src OPENBLAS_NUM_THREADS=1 .venv/bin/python src/report_trace_kappa_radial_audit.py`. The separate failed-endpoint probe is recorded explicitly in the JSON and is not part of the completed source-run settings.
