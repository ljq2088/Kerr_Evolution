# Schwarzschild alpha0.3 r20 finite-resolution run

First batch: scalar_batch_mg1_54daaacf7ded.json. Six scalar channels (2,2),(4,2),(6,2),
(0,0),(2,0),(4,0), from metric m1 and its conjugate. Parameters are
L18, nt18, nr8, h32 logarithmic horizon panel, source inner offset0.0005,
source outer320, Green inner offset0.0001, outer1000, two workers.

The massive |211> Schwarzschild cloud uses its complex eigenfrequency for the
radial eigenfunction and freezes temporal decay using its real frequency.
This explicitly approximate source is not an exact stationary Kerr threshold
cloud; the nonzero Klein-Gordon defect remains recorded in source metadata.

Boundary preflight schwarzschild_alpha03_r20_mg1_preflight.json checks the
six requested radial operators. Maximum last-term series ratio2.61203e-7
is a preflight diagnostic, not an outer-boundary convergence bound.

This is the first production-resolution batch for the Schwarzschild comparison,
not a completed Figure3 scan or a complete finite mode sum. Source, boundary,
angular and radial convergence remain required. Log:
outputs/schwarzschild_alpha03_r20_mg1.log.
