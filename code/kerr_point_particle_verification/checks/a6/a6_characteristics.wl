(* A6 targeted check: radial principal roots and boundary directions. *)
Module[
 {d, a, b, c, poly, roots, targetRoots, rootResiduals, rPlus, rH,
  a0, ah, scriCoeffResiduals, horizonCoeffResiduals, scriSpeeds,
  horizonSpeeds, scriSpeedResiduals, horizonSpeedResiduals,
  branchFormulaResiduals, angularPoly, angularRoots, angularRootResiduals,
  ass, horizonAss, discriminantResidual},
 d = 1 - yy^2;
 a = 8 mass (2 mass - spin^2 RR/len^2)
    (1 + 2 mass RR/len^2) - spin^2 d;
 b = -2 (len^2 - (8 mass^2 - spin^2) RR^2/len^2 +
    4 spin^2 mass RR^3/len^4);
 c = -(len^2 - 2 mass RR + spin^2 RR^2/len^2) RR^2/len^2;
 poly = a vv^2 - b vv + c;
 roots = vv /. Solve[poly == 0, vv];
 targetRoots = {b/2 + Sqrt[b^2/4 - a c],
    b/2 - Sqrt[b^2/4 - a c]}/a;
 rootResiduals = FullSimplify[poly /. vv -> #] & /@ targetRoots;
 discriminantResidual = FullSimplify[Discriminant[poly, vv] -
   (b^2 - 4 a c)];

 rPlus = mass + Sqrt[mass^2 - spin^2];
 rH = len^2/rPlus;
 ass = Element[{mass, spin, len, yy}, Reals] && mass > 0 && len > 0 &&
   -mass < spin < mass && -1 < yy < 1;
 horizonAss = ass && rPlus > 0;
 a0 = FullSimplify[a /. RR -> 0, ass];
 ah = FullSimplify[a /. RR -> rH, horizonAss];
 scriCoeffResiduals = FullSimplify[
   {c /. RR -> 0, (b /. RR -> 0) + 2 len^2}, ass];
 horizonCoeffResiduals = FullSimplify[
   {c /. RR -> rH, (b /. RR -> rH) - 4 mass rH}, horizonAss];
 scriSpeeds = {0, -2 len^2/a0};
 horizonSpeeds = {4 mass rH/ah, 0};
 scriSpeedResiduals = FullSimplify[
   {a0 scriSpeeds[[1]]^2 - (b /. RR -> 0) scriSpeeds[[1]] +
      (c /. RR -> 0),
    a0 scriSpeeds[[2]]^2 - (b /. RR -> 0) scriSpeeds[[2]] +
      (c /. RR -> 0)}, ass && a0 > 0];
 horizonSpeedResiduals = FullSimplify[
   {ah horizonSpeeds[[1]]^2 - (b /. RR -> rH) horizonSpeeds[[1]] +
      (c /. RR -> rH),
    ah horizonSpeeds[[2]]^2 - (b /. RR -> rH) horizonSpeeds[[2]] +
      (c /. RR -> rH)}, horizonAss && ah > 0];
 branchFormulaResiduals = {
   FullSimplify[{(bs/2 + Sqrt[bs^2/4])/as,
      (bs/2 - Sqrt[bs^2/4])/as} - {0, bs/as},
     Element[{as, bs}, Reals] && as > 0 && bs < 0],
   FullSimplify[{(bh/2 + Sqrt[bh^2/4])/ahh,
      (bh/2 - Sqrt[bh^2/4])/ahh} - {bh/ahh, 0},
     Element[{ahh, bh}, Reals] && ahh > 0 && bh > 0]};
 angularPoly = a vy^2 - d;
 angularRoots = {Sqrt[d/a], -Sqrt[d/a]};
 angularRootResiduals = FullSimplify[
   angularPoly /. vy -> #, ass && a > 0] & /@ angularRoots;
 <|"Stage" -> "A6 characteristics",
   "PrincipalPolynomial" -> "A v^2-B v+C",
   "RootFormulaResiduals" -> rootResiduals,
   "DiscriminantResidual" -> discriminantResidual,
   "SCRICoefficientResiduals" -> AssociationThread[
    {"C", "BPlus2L2"}, scriCoeffResiduals],
   "HorizonCoefficientResiduals" -> AssociationThread[
    {"C", "BMinus4MRH"}, horizonCoeffResiduals],
   "SCRISpeedResiduals" -> scriSpeedResiduals,
   "HorizonSpeedResiduals" -> horizonSpeedResiduals,
   "BoundaryBranchFormulaResiduals" -> branchFormulaResiduals,
   "DirectionChecks" -> <|
    "SCRITangent" -> TrueQ[scriSpeeds[[1]] == 0],
    "SCRINonzeroDecreasingR" -> FullSimplify[scriSpeeds[[2]] < 0,
      ass && a0 > 0],
    "HorizonNonzeroIncreasingR" -> FullSimplify[horizonSpeeds[[1]] > 0,
      horizonAss && ah > 0],
    "HorizonTangent" -> TrueQ[horizonSpeeds[[2]] == 0]|>,
   "AngularPrincipalPolynomial" -> "A vy^2-d",
   "AngularRoots" -> angularRoots,
   "AngularRootResiduals" -> angularRootResiduals,
   "AngularSpeedSquared" -> d/a,
   "CharacteristicsChecked" -> And @@ (# === 0 & /@ Flatten[{
     rootResiduals, discriminantResidual, scriCoeffResiduals,
     horizonCoeffResiduals, scriSpeedResiduals, horizonSpeedResiduals,
     branchFormulaResiduals, angularRootResiduals}]),
   "Assumptions" -> {"M>0,L>0,|a|<M,-1<y<1",
    "A(0,y)>0 and A(RH,y)>0 for radial branch ordering",
    "A(R,y)>0 and -1<y<1 for real angular roots"}|>
]
