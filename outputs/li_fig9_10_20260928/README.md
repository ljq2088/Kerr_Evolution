# Li v2 Figures 9 and 10: finite-resolution attempt

The run is active. Figure 9 has 20 of 24 visible paper points compared so far; four visible points need metric mg7--9. In addition, ell10,m2 is computed even though it is below the published plotting range. Figure 10 currently has every plotted channel for both clouds at r0=20M. Eleven orbital radii are queued: 20,18.3,4.1,10.1,30.1,41.1,42.1,50.1,18.1,18.5,40.1. This is a sparse scan, not a convergence-certified reconstruction of the entire continuous curves.

## Results at r0=20M

Figure9: 19/20 available plotted channels differ from the PDF marker readout by <=0.412%. The very weak ell8,m2 channel differs by -5.537% and requires precision/quadrature/outer-boundary checks. No normalization was fitted.

Figure10: the four cloud211 infinity modes differ from the published graph by a factor around 462. Multiplying the computed flux solely for diagnosis by the predetermined (r0*u^t)^2=469.035621857169 brings these four modes to within about 1.4--1.9%, but does not explain all horizon channels (ell2,m0 has about 51.5% residual). The reference itself is inconsistent in numerical scale across figures: Fig10/Fig9 for the same modes and r0=20 is 462.05--463.04. This is evidence of different displayed normalization, not proof of an author-code error.

For cloud322, the four infinity-mode reference/local ratios are approximately 1.03e8--1.05e8. Horizon discrepancies are not a single common factor. Applying (r0*u^t)^2 is insufficient. Figure10 is NOT reproduced.

## Independent controls

Both normalized cloud profiles match the original Fig3 vector curves at their vertices: cloud211 ratio median 1.000000014, cloud322 median 0.999995975. This rules out a simple overall cloud-amplitude error of 1e4 in the implemented background, relative to the published cloud profile.

For cloud322, separated and direct covariant source projections were compared for modes (4,4) and (1,1), at r=5.05175 and 24.12682. Relative differences are 2.5e-16--7.6e-15 for the SAME weighted spherical metric data. This tests the projection algebra, not the full correctness of the metric or Green boundary solution.

## Definitions and limitations

The physical effective orbital flux is used (including the cloud charge bookkeeping), in units epsilon^2 zeta^2. Unit particle-mass-squared/unit cloud-mass flux is divided by alpha^6. No r0*u^t factor is applied in any physical source, response, or flux output. The square-factor plots are separate diagnostics.

Input metric spheroidal L20 -> weighted spherical j18; angular metric q40, source q64, Fourier pmax12. Existing r20 data use nr12/h64 source quadrature; new orbits use nr24/h64 and phase-based panel refinement. Outer source radius320, source horizon offset5e-4; Green boundary order4. Complex spatial Leaver cloud with real temporal frequency, finite BL mass normalization. Authors' horizon prescription, higher precision and full convergence remain unverified.

Run and resume: PYTHONPATH=src OPENBLAS_NUM_THREADS=1 .venv/bin/python outputs/li_fig9_10_20260928/run_comparison.py --run --workers 3

The driver writes per-mode results, radial source checkpoints and versioned metric caches; refuses changed-code resumes; regenerates comparison plots after completed modes. Status is execution.json, details are group_r*_mg*.log. No time estimate is asserted for uncached metric groups.
