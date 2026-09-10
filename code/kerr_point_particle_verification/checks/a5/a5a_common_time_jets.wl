(* A5a: derive all common time jets from scalar definitions. *)
Module[
 {group, del, b, p, c, rGeo, vel, k, eta, fT, fPhi, uT, uR, uPhi,
  rDot, rDDot, rpC, rPdot, rPddot, phiDot, phiDDot,
  nHat, nExact, hatDeltaP, mBar, a0, dtr, dttR, dt, dtt,
  rPdotTarget, rPddotTarget,
  phiDDotTarget, nHatDotTarget, nHatDDotTarget, mDotTarget,
  mDDotTarget, lambda, lambdaDot, a0DotTarget, a0DDotTarget,
  amps, ampTargets, phaseTargets, xR, gR, gY, gauss, gamma1,
  gamma2, z, f, fTargets, residuals, ass},
 group = If[ValueQ[a5aGroup], a5aGroup, "Worldline"];
 del = rr^2 - 2 mass rr + spin^2;
 b = LZ - spin EE;
 p = EE (rr^2 + spin^2) - spin LZ;
 c = rr^2 + b^2;
 rGeo = p^2 - del c;
 vel = Sqrt[rGeo];
 k = c/(p + vel);
 eta = 2 + 4 mass/rr;
 fT = spin b + (rr^2 + spin^2) k + eta vel;
 fPhi = b + spin k;
 uT = fT/rr^2;
 uR = len^2 vel/rr^4;
 uPhi = fPhi/rr^2;
 rDot = -vel/(rr^2 uT);
 dtr[x_] := D[x, rr] rDot;
 dttR[x_] := D[dtr[x], rr] rDot;
 rDDot = dtr[rDot];
 rpC = len^2/rr;
 rPdot = dtr[rpC];
 rPddot = dttR[rpC];
 phiDot = uPhi/uT;
 phiDDot = dtr[phiDot];
 nHat = -k/(2 len^4);
 nExact = -del k/(2 rr^2);
 hatDeltaP = (len^2/rr)^2 del;
 mBar = -I b/(Sqrt[2] rr);
 a0 = mu len^2/(2 Pi rr^4 uT);

 (* Total T derivative on functions of radial position and orbital phase. *)
 dt[x_] := D[x, rr] rDot + D[x, ph] phiDot;
 dtt[x_] := dt[dt[x]];
 rPdotTarget = uR/uT;
 rPddotTarget = rDot (D[uR, rr] uT - uR D[uT, rr])/uT^2;
 phiDDotTarget = rDot (D[uPhi, rr] uT - uPhi D[uT, rr])/uT^2;
 nHatDotTarget = -D[k, rr] rDot/(2 len^4);
 nHatDDotTarget = -(D[k, {rr, 2}] rDot^2 + D[k, rr] rDDot)/
   (2 len^4);
 mDotTarget = I b rDot/(Sqrt[2] rr^2);
 mDDotTarget = -2 I b rDot^2/(Sqrt[2] rr^3) +
   I b rDDot/(Sqrt[2] rr^2);
 lambda = 4 rDot/rr + dtr[uT]/uT;
 lambdaDot = 4 (rDDot/rr - rDot^2/rr^2) +
   dttR[uT]/uT - (dtr[uT]/uT)^2;
 a0DotTarget = -lambda a0;
 a0DDotTarget = (lambda^2 - lambdaDot) a0;

 amps = {a0 nHat^2, a0 mBar nHat, a0 mBar^2};
 ampTargets = {
   {amps[[1]], a0DotTarget nHat^2 + 2 a0 nHat nHatDotTarget,
    a0DDotTarget nHat^2 + 4 a0DotTarget nHat nHatDotTarget +
     2 a0 (nHatDotTarget^2 + nHat nHatDDotTarget)},
   {amps[[2]], a0DotTarget mBar nHat +
     a0 (mDotTarget nHat + mBar nHatDotTarget),
    a0DDotTarget mBar nHat + 2 a0DotTarget
      (mDotTarget nHat + mBar nHatDotTarget) +
     a0 (mDDotTarget nHat + 2 mDotTarget nHatDotTarget +
       mBar nHatDDotTarget)},
   {amps[[3]], a0DotTarget mBar^2 + 2 a0 mBar mDotTarget,
    a0DDotTarget mBar^2 + 4 a0DotTarget mBar mDotTarget +
     2 a0 (mDotTarget^2 + mBar mDDotTarget)}};
 phaseTargets = Table[{
    ampTargets[[j, 1]] Exp[-I mmode ph],
    Exp[-I mmode ph] (ampTargets[[j, 2]] -
      I mmode phiDot ampTargets[[j, 1]]),
    Exp[-I mmode ph] (ampTargets[[j, 3]] -
      2 I mmode phiDot ampTargets[[j, 2]] -
      I mmode phiDDot ampTargets[[j, 1]] -
      mmode^2 phiDot^2 ampTargets[[j, 1]])}, {j, 3}];

 xR = rad - rpC;
 gR = Exp[-xR^2/(2 sigR^2)]/(Sqrt[2 Pi] sigR);
 gY = Exp[-yy^2/(2 sigY^2)]/(Sqrt[2 Pi] sigY);
 gauss = gR gY;
 gamma1 = xR rPdot/sigR^2;
 gamma2 = gamma1^2 + (xR rPddot - rPdot^2)/sigR^2;
 z = zAmp[rr] Exp[-I mmode ph];
 f = z gauss;
 fTargets = {f,
   (dt[z] + z gamma1) gauss,
   (dtt[z] + 2 dt[z] gamma1 + z gamma2) gauss};

 ass = Element[{mass, spin, len, mu, EE, LZ, rr, ph, mmode,
    sigR, sigY, rad, yy}, Reals] && mass > 0 && len > 0 && rr > 0 &&
   sigR > 0 && sigY > 0 && rGeo > 0 && p + vel != 0 && uT > 0;
 residuals = Which[group === "Worldline", <|
    "Rdot" -> FullSimplify[rPdot - rPdotTarget, ass],
    "Rddot" -> FullSimplify[rPddot - rPddotTarget, ass],
    "Phiddot" -> FullSimplify[phiDDot - phiDDotTarget, ass],
    "NhatDot" -> FullSimplify[dtr[nHat] - nHatDotTarget, ass],
    "NhatDDot" -> FullSimplify[dttR[nHat] - nHatDDotTarget, ass],
    "MbarDot" -> FullSimplify[dtr[mBar] - mDotTarget, ass],
    "MbarDDot" -> FullSimplify[dttR[mBar] - mDDotTarget, ass],
    "A0Dot" -> FullSimplify[dtr[a0] - a0DotTarget, ass],
    "A0DDot" -> FullSimplify[dttR[a0] - a0DDotTarget, ass],
    "ExtensionRecovery" -> FullSimplify[{
      hatDeltaP nHat - nExact,
      hatDeltaP^2 a0 nHat^2 - a0 nExact^2,
      hatDeltaP a0 mBar nHat - a0 mBar nExact,
      a0 mBar^2 - a0 mBar^2}, ass]|>,
   group === "AmplitudePhase", <|
    "AmplitudeJets" -> Table[FullSimplify[{
      amps[[j]] - ampTargets[[j, 1]],
      dtr[amps[[j]]] - ampTargets[[j, 2]],
      dttR[amps[[j]]] - ampTargets[[j, 3]]}, ass], {j, 3}],
    "PhaseJets" -> Table[FullSimplify[{
      amps[[j]] Exp[-I mmode ph] - phaseTargets[[j, 1]],
      dt[amps[[j]] Exp[-I mmode ph]] - phaseTargets[[j, 2]],
      dtt[amps[[j]] Exp[-I mmode ph]] - phaseTargets[[j, 3]]}, ass],
      {j, 3}]|>,
   group === "Gaussian", <|
    "FullLineNormalization" -> FullSimplify[{
      Integrate[Exp[-(x - rpC)^2/(2 sigR^2)]/(Sqrt[2 Pi] sigR),
        {x, -Infinity, Infinity}, Assumptions -> sigR > 0] - 1,
      Integrate[Exp[-x^2/(2 sigY^2)]/(Sqrt[2 Pi] sigY),
        {x, -Infinity, Infinity}, Assumptions -> sigY > 0] - 1}, ass],
    "CenteredFirstMoments" -> FullSimplify[{
      Integrate[(x - rpC) Exp[-(x - rpC)^2/(2 sigR^2)]/
        (Sqrt[2 Pi] sigR), {x, -Infinity, Infinity},
        Assumptions -> sigR > 0],
      Integrate[x Exp[-x^2/(2 sigY^2)]/(Sqrt[2 Pi] sigY),
        {x, -Infinity, Infinity}, Assumptions -> sigY > 0]}, ass],
    "GaussianJets" -> FullSimplify[{
      gauss - fTargets[[1]]/z,
      dtr[gauss] - gamma1 gauss,
      dttR[gauss] - gamma2 gauss}, ass],
    "ProductGaussianJets" -> FullSimplify[{
      f - fTargets[[1]], dt[f] - fTargets[[2]],
      dtt[f] - fTargets[[3]]}, ass]|>, True, <||>];
 <|"Stage" -> "A5a common time jets",
   "ExecutionGroup" -> group,
   "Residuals" -> residuals,
   "AllRequiredResidualsZero" ->
    And @@ (# === 0 & /@ Flatten[Values[residuals]]),
   "UndefinedDotQuantitiesUsedAsInputs" -> False,
   "ScalarInputs" ->
    {"M", "a", "L", "mu", "m", "E", "Lz", "rp", "PhiP",
     "sigmaR", "sigmaY", "R", "y"},
   "Assumptions" -> {
    "fixed positive widths", "Rgeo>0 inward branch",
    "P+Sqrt[Rgeo]!=0,uT>0", "equatorial yp=0"},
   "Excluded" -> {"turning point", "time-dependent widths",
    "finite-domain normalization", "axis completion"}|>
]
