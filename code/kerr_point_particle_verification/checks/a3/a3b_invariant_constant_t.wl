(* A3b route 1: invariant stress tensor reduced on a constant-T slice. *)
Module[
 {del, sig, gBL, jac, inverseDet, metricDetRes, inverseDetRes,
  sqrtgHDerived, sqrtgHTarget, volumeRes, tauRoot, rootConditionRes,
  rootDerivative, rootDerivativeAtRoot, absBranchRes, tauJacobian,
  tauJacobianRes, invariantWeight, targetWeight, weightRes, uVec,
  dist3, tensorDerived, tensorTarget, tensorResiduals, required, ass},
 del = rr^2 - 2 mass rr + spin^2;
 sig = rr^2 + spin^2 Cos[th]^2;
 gBL = {{-(1 - 2 mass rr/sig), 0, 0,
     -2 mass spin rr Sin[th]^2/sig},
   {0, sig/del, 0, 0}, {0, 0, sig, 0},
   {-2 mass spin rr Sin[th]^2/sig, 0, 0,
    Sin[th]^2 ((rr^2 + spin^2)^2 -
      spin^2 del Sin[th]^2)/sig}};
 jac = {{1, hP, 0, 0}, {0, -ll^2/rr^2, 0, 0},
   {0, 0, 1, 0}, {0, qP, 0, 1}};
 ass = Element[{mass, spin, rr, ll, th, uT}, Reals] &&
   mass > 0 && ll > 0 && rr > 0 && 0 < th < Pi && uT > 0 &&
   del != 0;
 metricDetRes = FullSimplify[Det[gBL] + sig^2 Sin[th]^2, ass];
 inverseDet = FullSimplify[Abs[1/Det[jac]], ass];
 inverseDetRes = FullSimplify[inverseDet - rr^2/ll^2, ass];
 sqrtgHDerived = FullSimplify[Sqrt[-Det[gBL]] inverseDet, ass];
 sqrtgHTarget = sig rr^2 Sin[th]/ll^2;
 volumeRes = FullSimplify[sqrtgHDerived - sqrtgHTarget, ass];

 tauRoot = Tobs - tWorld[tau];
 rootConditionRes = tauRoot /. tau -> tauP /. tWorld[tauP] -> Tobs;
 rootDerivative = D[tauRoot, tau];
 rootDerivativeAtRoot = rootDerivative /. tau -> tauP /.
   Derivative[1][tWorld][tauP] -> uT;
 absBranchRes = FullSimplify[Abs[rootDerivativeAtRoot] - uT, ass];
 tauJacobian = FullSimplify[1/Abs[rootDerivativeAtRoot], ass];
 tauJacobianRes = FullSimplify[tauJacobian - 1/uT, ass];

 invariantWeight = FullSimplify[mu tauJacobian/sqrtgHDerived, ass];
 targetWeight = mu/(sqrtgHTarget uT);
 weightRes = FullSimplify[invariantWeight - targetWeight, ass];
 uVec = {uT, uR, uTheta, uPhi};
 dist3 = deltaR deltaTheta deltaPhi;
 tensorDerived = invariantWeight Outer[Times, uVec, uVec] dist3;
 tensorTarget = targetWeight Outer[Times, uVec, uVec] dist3;
 tensorResiduals = FullSimplify[tensorDerived - tensorTarget, ass];
 required = Flatten[{metricDetRes, inverseDetRes, volumeRes,
   rootConditionRes, absBranchRes, tauJacobianRes, weightRes,
   tensorResiduals}];

 <|"Stage" -> "A3b invariant constant-T route",
   "PhysicalStart" ->
    "mu Integral[u^A u^B delta^4(X-z(tau))/sqrt(-g_H),{tau}]",
   "GeometryResiduals" -> <|
    "BLMetricDeterminant" -> metricDetRes,
    "InverseCoordinateJacobianDeterminant" -> inverseDetRes,
    "HyperboloidalVolume" -> volumeRes|>,
   "TauRootEvidence" -> <|
    "RootFunction" -> tauRoot,
    "RootConditionResidual" -> rootConditionRes,
    "RootDerivative" -> rootDerivative,
    "AbsoluteValueBranchResidual" -> absBranchRes,
    "OneOverAbsDerivativeResidual" -> tauJacobianRes|>,
   "RemainingSpatialDeltas" ->
    {"delta(R-Rp(T))", "delta(theta-thetaP(T))", "delta(Phi-PhiP(T))"},
   "DistributionWeightResidual" -> weightRes,
   "StressTensorEndpointResiduals" -> tensorResiduals,
   "Assumptions" -> {
    "C1 worldline near the root",
    "one simple transverse local root with uT!=0",
    "future branch uT>0",
    "M>0,L>0,r>0,0<theta<Pi and nonsingular coordinate Jacobian"},
   "GlobalWorldlineMonotonicityVerified" -> False,
   "GlobalUniqueIntersectionVerified" -> False,
   "LocalSimpleTransverseRootAssumed" -> True,
   "AllRequiredResidualsZero" -> And @@ (# === 0 & /@ required)|>
]
