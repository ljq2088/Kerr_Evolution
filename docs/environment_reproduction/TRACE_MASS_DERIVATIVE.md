# Differentiated forced trace resolvent and kappa jets

Diagnostic implementation: `environment_trace_variation.py`. Production imports
are unchanged, so ongoing metric caches and source calculations are unaffected.

For a particle at r0, the separated trace source is J=-16 pi S(pi/2)/ut.
Writing nu for auxiliary mass squared, J_nu=-16 pi S_nu(pi/2)/ut.
The retarded coefficients are

    z_infinity = J R_in(r0)/W,
    z_horizon  = J R_up(r0)/W,
    z_nu = (J_nu R + J R_nu)/W - z W_nu/W.

Inside the orbit F=z_horizon R_in and outside F=z_infinity R_up.
The field derivative is h_nu=F_nu S+F S_nu. The implementation includes
all these terms; holding the source projector fixed would be a different
finite-multipole resolvent. Continuity and radial derivative jumps are

    [F]=[F_nu]=0,
    Delta0 [F_r]=J,   Delta0 [F_rnu]=J_nu.

The local radial Taylor recurrence is differentiated using
Delta D''+Delta' D'+V D=(r^2+A_nu)R. The angular recurrence uses

    E''+cot(theta) E'+U E=(a^2 cos(theta)^2-A_nu)S,
    U=a^2(omega^2-nu)cos(theta)^2+A-m^2/sin(theta)^2.

Thus both the original trace and its derivative are available as full
mixed Taylor jets. At nu=0, kappa=-i omega h_nu obeys Box kappa=-i omega h
away from the particle. This follows by differentiating (Box-nu)h_nuparam=0;
the source distribution and angular sum require the jump/projector treatment.

Verification uses the independent existing finite-difference resolvent at
three orbits/modes, including negative m, on both sides of the particle.
Relative discrepancies are below 2.1e-12 for fields and first derivatives,
and below 1.2e-8 for their mass derivatives. Particle jumps pass separately.
Eighth-order jets at r=4.5 and 12 for r0=10, a=.877, ell=3, m=1 agree with
the existing kappa jets to 4.5e-9 in the maximum-coefficient norm. The covariant
wave equation residual relative to its target is below 1.4e-14 at those points.
These finite checks do not certify high-L metric reconstruction or reproduce
paper fluxes. The next check is substituting this independent kappa into the
metric source-matching diagnostic, including its high-L sequence, before
considering any production change. Outer boundary convergence remains separate.
