# Li field assembler independent review — 2026-09-17

Reviewed `src/report_li_field_alignment_plot.py` against the actual calibrated reference JSON/NPZ and the current `EnvironmentalWake`, `SampledResponse`, and provenance validator. No radial/source calculation was rerun. Production modules were not edited.

## Confirmed

- NPZ keys, array order `(orbit,radius,phi)=(2,9,72)`, and orbit lookup match the generated reference. Every cell is finite and RGB distance is zero; intervals contain their nominal midpoint.
- Full mode set has 18 reflection-allowed modes for scalar ell=2..5; the separate six positive modes are `(2,2),(3,3),(4,2),(4,4),(5,3),(5,5)`. They preserve the paper's caption/body ambiguity rather than choosing by fit.
- `EnvironmentalWake` produces the physical scalar field with the cached cloud mass. Dividing by `alpha^3 sqrt(cloud_mass)` gives the specified dimensionless expansion coefficient. All current inputs have mass=1, so explicitly including sqrt(mass) does not change their numbers.
- Local phi=0 is right, phi=pi/2 is up, and time defaults to zero. `imshow(origin=upper,extent=(-200,200,-200,200))` keeps the PNG top at positive y. The reference's conventional azimuth registration remains an assumption documented by the digitizer; no phase or rotation is fitted.
- Ring color intervals are applied pointwise. Angular RMS interval bounds follow monotonicity for nonnegative amplitudes. These discrete-azimuth RMS values are not claimed to reconstruct missing complex phase.
- Source samples are checked against the current implementation and their sample checksum; `SampledResponse.from_report` separately requires a completed finite response. Mixing numerical discretizations is explicit.

## Report-only fixes

The original plotting code used continuous standard `viridis` for local fields, but the published colorbar has a discrete, nonuniform mapping. At its midpoint the original RGB is `[50,182,122]`, whereas standard viridis gives `[33,144,141]`. That would make equal numerical amplitudes look different. The plot now uses the **exact original colorbar RGB rows**, reversed into ascending field order, with the same nominal upper endpoint .37. This changes display colors only; it does not rescale field values.

Also added:

- explicit complete-response, expected alpha/spin, and positive cloud amplitude checks;
- SHA matching of the actual PNG and NPZ to the digitization JSON, plus NPZ shape/finite/color-interval validation;
- the missing general sqrt(cloud_mass) in all three field evaluations (neutral for this mass=1 run);
- point-table lower/upper quantization intervals and inclusion flag, in JSON and CSV;
- source/assembly dependency hash consistency at the end of evaluation, with dependency hashes recorded in the output.

## Validation scope

Syntax parsed successfully. Actual reference SHA, all NPZ cells and intervals, and construction of the 864-entry original LUT passed cheap checks. Numerical wake assembly remains for the parent process after all fresh source modes are complete; this review does not claim to have run it.
