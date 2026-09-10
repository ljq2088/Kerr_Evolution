(* A3c: Kinnersley covectors, contractions, and the three A3 projections. *)
Module[
 {del, sig, b, p, c, rGeo, sqrtR, k, hPrime, qPrime, gBL,
  nUp, mbUp, nLow, mbLow, jac, invJac, nH, mbH, nHTarget,
  mbHTarget, uBL, uH, ut, ur, uphi, uT, uR, uPhi, geoAss,
  nLowerRes, mbLowerRes, nPullRes, mbPullRes, nBLContract,
  mbBLContract, nHContract, mbHContract, nTarget, mbTarget,
  blContractionRes, hContractionRes, contractionRouteRes,
  sqrtgHEq, weightInvariant, weightBL, weightEndpoint,
  invariantWeightRes, blWeightRes, dist3, projInvariant, projBL,
  projEndpoint, invariantProjectionRes, blProjectionRes,
  routeProjectionRes, horizonNScalingRes, factorizationDiagnostic,
  required},
 del = rr^2 - 2 mass rr + spin^2;
 sig = rr^2 + spin^2 Cos[th]^2;
 b = LZ - spin EE;
 p = EE (rr^2 + spin^2) - spin LZ;
 c = rr^2 + b^2;
 rGeo = p^2 - del c;
 sqrtR = Sqrt[rGeo];
 k = c/(p + sqrtR);
 hPrime = (rr^2 + spin^2)/del - 2 - 4 mass/rr;
 qPrime = spin/del;
 gBL = {{-(1 - 2 mass rr/sig), 0, 0,
     -2 mass spin rr Sin[th]^2/sig},
   {0, sig/del, 0, 0}, {0, 0, sig, 0},
   {-2 mass spin rr Sin[th]^2/sig, 0, 0,
    Sin[th]^2 ((rr^2 + spin^2)^2 -
      spin^2 del Sin[th]^2)/sig}};
 nUp = {rr^2 + spin^2, -del, 0, spin}/(2 sig);
 mbUp = {-I spin Sin[th], 0, 1, -I Csc[th]}/
   (Sqrt[2] (rr - I spin Cos[th]));
 geoAss = Element[{mass, spin, rr, EE, LZ, ll}, Reals] &&
   mass > 0 && -mass < spin < mass && ll > 0 && rr > 0 &&
   del > 0 && rGeo >= 0 && p + sqrtR != 0;
 nLow = FullSimplify[gBL.nUp, geoAss && 0 < th < Pi];
 mbLow = FullSimplify[gBL.mbUp, geoAss && 0 < th < Pi];
 jac = {{1, hPrime, 0, 0}, {0, -ll^2/rr^2, 0, 0},
   {0, 0, 1, 0}, {0, qPrime, 0, 1}};
 invJac = FullSimplify[Inverse[jac], geoAss];
 nH = FullSimplify[Transpose[invJac].nLow, geoAss && 0 < th < Pi];
 mbH = FullSimplify[Transpose[invJac].mbLow, geoAss && 0 < th < Pi];
 nHTarget = {-del/(2 sig),
   del rr^2/(ll^2 sig) (1 + 2 mass/rr), 0,
   spin del Sin[th]^2/(2 sig)};
 mbHTarget = {I spin Sin[th],
   -I spin rr^2 Sin[th]/ll^2 (2 + 4 mass/rr), sig,
   -I (rr^2 + spin^2) Sin[th]}/
   (Sqrt[2] (rr - I spin Cos[th]));
 nLowerRes = FullSimplify[nLow -
   {-del/(2 sig), -1/2, 0, spin del Sin[th]^2/(2 sig)},
   geoAss && 0 < th < Pi];
 mbLowerRes = FullSimplify[mbLow -
   {I spin Sin[th], 0, sig, -I (rr^2 + spin^2) Sin[th]}/
    (Sqrt[2] (rr - I spin Cos[th])), geoAss && 0 < th < Pi];
 nPullRes = FullSimplify[nH - nHTarget, geoAss && 0 < th < Pi];
 mbPullRes = FullSimplify[mbH - mbHTarget, geoAss && 0 < th < Pi];

 ut = (spin b + (rr^2 + spin^2) p/del)/rr^2;
 ur = -sqrtR/rr^2;
 uphi = (b + spin p/del)/rr^2;
 uT = (spin b + (rr^2 + spin^2) k +
    (2 + 4 mass/rr) sqrtR)/rr^2;
 uR = ll^2 sqrtR/rr^4;
 uPhi = (b + spin k)/rr^2;
 uBL = {ut, ur, 0, uphi};
 uH = {uT, uR, 0, uPhi};
 nBLContract = FullSimplify[(nLow /. th -> Pi/2).uBL, geoAss];
 mbBLContract = FullSimplify[(mbLow /. th -> Pi/2).uBL, geoAss];
 nHContract = FullSimplify[(nH /. th -> Pi/2).uH, geoAss];
 mbHContract = FullSimplify[(mbH /. th -> Pi/2).uH, geoAss];
 nTarget = -del k/(2 rr^2);
 mbTarget = -I b/(Sqrt[2] rr);
 blContractionRes = FullSimplify[
   {nBLContract - nTarget, mbBLContract - mbTarget}, geoAss];
 hContractionRes = FullSimplify[
   {nHContract - nTarget, mbHContract - mbTarget}, geoAss];
 contractionRouteRes = FullSimplify[
   {nBLContract - nHContract, mbBLContract - mbHContract}, geoAss];

 sqrtgHEq = rr^4/ll^2;
 weightInvariant = mu/(sqrtgHEq uT);
 weightBL = FullSimplify[
   mu/(rr^2 ut) (ut/uT) (ll^2/rr^2), geoAss && ut > 0 && uT > 0];
 weightEndpoint = mu ll^2/(rr^4 uT);
 invariantWeightRes = FullSimplify[weightInvariant - weightEndpoint,
   geoAss && uT > 0];
 blWeightRes = FullSimplify[weightBL - weightEndpoint,
   geoAss && ut > 0 && uT > 0];
 dist3 = deltaR deltaTheta deltaPhi;
 projInvariant = weightInvariant
   {nHContract^2, mbHContract nHContract, mbHContract^2} dist3;
 projBL = weightBL
   {nBLContract^2, mbBLContract nBLContract, mbBLContract^2} dist3;
 projEndpoint = weightEndpoint
   {nTarget^2, mbTarget nTarget, mbTarget^2} dist3;
 invariantProjectionRes = FullSimplify[
   projInvariant - projEndpoint, geoAss && uT > 0];
 blProjectionRes = FullSimplify[
   projBL - projEndpoint, geoAss && ut > 0 && uT > 0];
 routeProjectionRes = FullSimplify[
   projInvariant - projBL, geoAss && ut > 0 && uT > 0];
 horizonNScalingRes = FullSimplify[nTarget/del + k/(2 rr^2), geoAss];
 factorizationDiagnostic = FullSimplify[
   {nHContract^2 - nTarget^2,
    mbHContract nHContract - mbTarget nTarget,
    mbHContract^2 - mbTarget^2}, geoAss];
 required = Flatten[{nLowerRes, mbLowerRes, nPullRes, mbPullRes,
   blContractionRes, hContractionRes, contractionRouteRes,
   invariantWeightRes, blWeightRes, invariantProjectionRes,
   blProjectionRes, routeProjectionRes, horizonNScalingRes}];

 <|"Stage" -> "A3c tetrad projections",
   "PrimaryStart" ->
    "P08 Kinnersley vectors plus BL metric, coordinate pullback, and the two verified T_H weights",
   "CovectorResiduals" -> <|
    "BLn" -> nLowerRes, "BLmbar" -> mbLowerRes,
    "HyperboloidalN" -> nPullRes,
    "HyperboloidalMbar" -> mbPullRes|>,
   "ContractionEndpointResiduals" -> <|
    "BLRoute" -> AssociationThread[{"Np", "Mbarp"}, blContractionRes],
    "InvariantRoute" -> AssociationThread[{"Np", "Mbarp"}, hContractionRes],
    "RouteComparison" -> AssociationThread[{"Np", "Mbarp"},
      contractionRouteRes]|>,
   "TheoryContractions" -> <|
    "Np" -> "-Delta K/(2 r^2)",
    "Mbarp" -> "-I (Lz-aE)/(Sqrt[2] r)"|>,
   "RouteWeightResiduals" -> <|
    "Invariant" -> invariantWeightRes, "BLPullback" -> blWeightRes|>,
   "ProjectionEndpointResiduals" -> <|
    "InvariantRoute" -> AssociationThread[
      {"Tnn", "TmbarN", "TmbarMbar"}, invariantProjectionRes],
    "BLPullbackRoute" -> AssociationThread[
      {"Tnn", "TmbarN", "TmbarMbar"}, blProjectionRes],
    "WholeRouteComparison" -> AssociationThread[
      {"Tnn", "TmbarN", "TmbarMbar"}, routeProjectionRes]|>,
   "NpOverDeltaRegularScalingResidual" -> horizonNScalingRes,
   "Diagnostics" -> <|
    "RankOneProjectionFactorizationResiduals" -> factorizationDiagnostic,
    "Classification" -> "same-definition diagnostic, not standalone pass evidence"|>,
   "Assumptions" -> {
    "equatorial Q=0 inward geodesic, Rgeo>=0",
    "M>0,|a|<M,L>0,r>0,Delta>0,P+Sqrt[Rgeo]!=0",
    "ut>0,uT>0 for the two distribution weights",
    "three common deltas are relative to dR dtheta dPhi"},
   "A2FourBlockOperatorsUsed" -> False,
   "AllRequiredResidualsZero" -> And @@ (# === 0 & /@ required)|>
]
