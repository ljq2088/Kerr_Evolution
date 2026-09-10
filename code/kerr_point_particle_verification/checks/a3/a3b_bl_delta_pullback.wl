(* A3b route 2: pull P08 Eq. 26 from BL to (T,R,theta,Phi). *)
Module[
 {jac, uBL, uH, root, rootConditionRes, rootDerivative,
  rootDerivativeAtSupport, ass, rootAbsRes, rootJacobian,
  rootJacobianRes, radialDerivative, radialJacobian, radialJacobianRes,
  thetaJacobianRes, azArgument, azSupportRes, azDerivative, azAbsRes,
  blWeightPulled, targetWeight, weightRes, dist3, tensorPulled,
  tensorTarget, tensorResiduals, velocityResiduals, required},
 jac = {{1, hp, 0, 0}, {0, -ll^2/rp^2, 0, 0},
   {0, 0, 1, 0}, {0, qp, 0, 1}};
 uBL = {ut, ur, uTheta, uPhiBL};
 uH = jac.uBL;
 ass = Element[{rp, ll, thetaP, ut, ur, hp, uT, sigP}, Reals] &&
   rp > 0 && ll > 0 && 0 < thetaP < Pi && sigP > 0 &&
   ut > 0 && uT > 0 && uT == ut + hp ur;

 root = rVar - rWorld[Tobs - hWorld[rVar]];
 rootConditionRes = root /. rVar -> rp /.
   Tobs - hWorld[rp] -> tP /. rWorld[tP] -> rp;
 rootDerivative = D[root, rVar];
 rootDerivativeAtSupport = rootDerivative /. rVar -> rp /.
   Tobs - hWorld[rp] -> tP /.
   Derivative[1][rWorld][tP] -> ur/ut /.
   Derivative[1][hWorld][rp] -> hp;
 rootAbsRes = FullSimplify[
   Abs[rootDerivativeAtSupport] - uT/ut, ass];
 rootJacobian = FullSimplify[1/Abs[rootDerivativeAtSupport], ass];
 rootJacobianRes = FullSimplify[rootJacobian - ut/uT, ass];
 radialDerivative = D[ll^2/rVar, rVar] /. rVar -> rp;
 radialJacobian = FullSimplify[Abs[radialDerivative], ass];
 radialJacobianRes = FullSimplify[radialJacobian - ll^2/rp^2, ass];
 thetaJacobianRes = 1 - 1;

 azArgument = PhiObs - qWorld[rVar] -
   phiWorld[Tobs - hWorld[rVar]];
 azSupportRes = azArgument /. rVar -> rp /.
   Tobs - hWorld[rp] -> tP /.
   phiWorld[tP] -> PhiP - qWorld[rp];
 azSupportRes = FullSimplify[azSupportRes - (PhiObs - PhiP)];
 azDerivative = D[azArgument, PhiObs];
 azAbsRes = FullSimplify[Abs[azDerivative] - 1];

 (* Eq. 26 weight times all spatial delta Jacobians. *)
 blWeightPulled = FullSimplify[
   mu/(sigP Sin[thetaP] ut) rootJacobian radialJacobian, ass];
 targetWeight = mu ll^2/(sigP rp^2 Sin[thetaP] uT);
 weightRes = FullSimplify[blWeightPulled - targetWeight, ass];
 velocityResiduals = FullSimplify[uH -
   {ut + hp ur, -ll^2 ur/rp^2, uTheta, uPhiBL + qp ur}, ass];
 dist3 = deltaR deltaTheta deltaPhi;
 tensorPulled = blWeightPulled Outer[Times, uH, uH] dist3;
 tensorTarget = targetWeight Outer[Times,
   {uT, -ll^2 ur/rp^2, uTheta, uPhiBL + qp ur},
   {uT, -ll^2 ur/rp^2, uTheta, uPhiBL + qp ur}] dist3;
 tensorResiduals = FullSimplify[tensorPulled - tensorTarget, ass];
 required = Flatten[{rootConditionRes, rootAbsRes, rootJacobianRes,
   radialJacobianRes, thetaJacobianRes, azSupportRes, azAbsRes,
   weightRes, velocityResiduals, tensorResiduals}];

 <|"Stage" -> "A3b P08 Eq. 26 BL delta pullback",
   "PhysicalStart" ->
    "mu u^mu u^nu delta(r-rp(t)) delta(theta-thetaP) delta(phi-phiP(t))/(Sigma sin(theta) ut)",
   "RadialRootEvidence" -> <|
    "RootFunction" -> root,
    "RootConditionResidual" -> rootConditionRes,
    "RootDerivative" -> rootDerivative,
    "AbsoluteValueResidual" -> rootAbsRes,
    "OneOverAbsDerivativeResidual" -> rootJacobianRes,
    "RCoordinateJacobianResidual" -> radialJacobianRes|>,
   "AngularDeltaJacobianResidual" -> thetaJacobianRes,
   "AzimuthalDeltaEvidence" -> <|
    "Argument" -> azArgument,
    "SupportResidual" -> azSupportRes,
    "DerivativeWithRespectToPhi" -> azDerivative,
    "AbsoluteValueResidual" -> azAbsRes|>,
   "VelocityJacobianResiduals" -> velocityResiduals,
   "DistributionWeightResidual" -> weightRes,
   "StressTensorEndpointResiduals" -> tensorResiduals,
   "EndpointSpatialDeltas" ->
    {"delta(R-Rp(T))", "delta(theta-thetaP)", "delta(Phi-PhiP(T))"},
   "Assumptions" -> {
    "C1 radial and azimuthal worldline functions near support",
    "one simple local radial root",
    "ut>0,uT>0,uT=ut+hp ur",
    "rp>0,L>0,0<thetaP<Pi,SigmaP>0",
    "azimuthal equality is local modulo 2Pi"},
   "GlobalRadialRootUniquenessVerified" -> False,
   "AzimuthalGlobalBranchVerified" -> False,
   "LocalRootAndSupportComparison" -> True,
   "AllRequiredResidualsZero" -> And @@ (# === 0 & /@ required)|>
]
