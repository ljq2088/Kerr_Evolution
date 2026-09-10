(* A3a: equatorial Kerr geodesic and horizon-regular hyperboloidal velocity. *)
Module[
 {del, rPlus, b, p, c, rGeo, sqrtR, gEq, solution, ut, uphi,
  urSquared, ur, utTarget, uphiTarget, geoAss, energyRes, lzRes,
  normRes, urSquaredRes, inwardRes, hPrime, qPrime, kSingular,
  kRegular, kRes, uTRaw, uRRaw, uPhiRaw, uTTarget, uRTarget,
  uPhiTarget, uTRes, uRRes, uPhiRes, pPlus, cPlus, horizonBaseAss,
  horizonAss, deltaPlusRes, rGeoPlusRes, pPlusDefinitionRes,
  sqrtBranchRes, kPlusRes,
  required},
 del = rr^2 - 2 mass rr + spin^2;
 rPlus = mass + Sqrt[mass^2 - spin^2];
 b = LZ - spin EE;
 p = EE (rr^2 + spin^2) - spin LZ;
 c = rr^2 + b^2;
 rGeo = p^2 - del c;
 sqrtR = Sqrt[rGeo];
 gEq = {{-(1 - 2 mass/rr), 0, 0, -2 mass spin/rr},
   {0, rr^2/del, 0, 0}, {0, 0, rr^2, 0},
   {-2 mass spin/rr, 0, 0,
    rr^2 + spin^2 + 2 mass spin^2/rr}};

 solution = FullSimplify[LinearSolve[
   {{-gEq[[1, 1]], -gEq[[1, 4]]},
    {gEq[[4, 1]], gEq[[4, 4]]}}, {EE, LZ}], del != 0];
 ut = solution[[1]];
 uphi = solution[[2]];
 urSquared = FullSimplify[
   (-1 - {ut, 0, 0, uphi}.gEq.{ut, 0, 0, uphi})/gEq[[2, 2]]];
 ur = -sqrtR/rr^2;
 utTarget = (spin b + (rr^2 + spin^2) p/del)/rr^2;
 uphiTarget = (b + spin p/del)/rr^2;
 geoAss = Element[{mass, spin, rr, EE, LZ, len}, Reals] &&
   mass > 0 && -mass < spin < mass && rr > rPlus && len > 0 &&
   del > 0 && rGeo >= 0;
 energyRes = FullSimplify[-gEq[[1]].{ut, ur, 0, uphi} - EE, geoAss];
 lzRes = FullSimplify[gEq[[4]].{ut, ur, 0, uphi} - LZ, geoAss];
 normRes = FullSimplify[{ut, ur, 0, uphi}.gEq.{ut, ur, 0, uphi} + 1,
   geoAss];
 urSquaredRes = FullSimplify[urSquared - rGeo/rr^4, geoAss];
 inwardRes = FullSimplify[ur + sqrtR/rr^2, geoAss];

 hPrime = (rr^2 + spin^2)/del - 2 - 4 mass/rr;
 qPrime = spin/del;
 kSingular = (p - sqrtR)/del;
 kRegular = c/(p + sqrtR);
 kRes = FullSimplify[Together[kSingular - kRegular],
   geoAss && p + sqrtR != 0];
 uTRaw = FullSimplify[ut + hPrime ur, geoAss];
 uRRaw = FullSimplify[-len^2 ur/rr^2, geoAss];
 uPhiRaw = FullSimplify[uphi + qPrime ur, geoAss];
 uTTarget = (spin b + (rr^2 + spin^2) kRegular +
    (2 + 4 mass/rr) sqrtR)/rr^2;
 uRTarget = len^2 sqrtR/rr^4;
 uPhiTarget = (b + spin kRegular)/rr^2;
 uTRes = FullSimplify[Together[uTRaw - uTTarget],
   geoAss && p + sqrtR != 0];
 uRRes = FullSimplify[uRRaw - uRTarget, geoAss];
 uPhiRes = FullSimplify[Together[uPhiRaw - uPhiTarget],
   geoAss && p + sqrtR != 0];

 pPlus = FullSimplify[p /. rr -> rPlus,
   Element[{mass, spin, EE, LZ}, Reals] && mass > 0 && -mass < spin < mass];
 cPlus = c /. rr -> rPlus;
 horizonBaseAss = Element[{mass, spin, EE, LZ}, Reals] &&
   mass > 0 && -mass < spin < mass;
 horizonAss = Element[pH, Reals] && pH > 0;
 deltaPlusRes = FullSimplify[del /. rr -> rPlus, horizonBaseAss];
 rGeoPlusRes = FullSimplify[(rGeo /. rr -> rPlus) - pPlus^2,
   horizonBaseAss];
 pPlusDefinitionRes = FullSimplify[
   pPlus - (2 mass EE rPlus - spin LZ), horizonBaseAss];
 sqrtBranchRes = FullSimplify[Sqrt[pH^2] - pH, horizonAss];
 kPlusRes = FullSimplify[
   cPlus/(pH + Sqrt[pH^2]) - cPlus/(2 pH), horizonAss];

 required = {FullSimplify[ut - utTarget, geoAss],
   FullSimplify[uphi - uphiTarget, geoAss], urSquaredRes, inwardRes,
   energyRes, lzRes, normRes, kRes, uTRes, uRRes, uPhiRes,
   deltaPlusRes, rGeoPlusRes, pPlusDefinitionRes, sqrtBranchRes, kPlusRes};
 <|"Stage" -> "A3a geodesic velocities",
   "PrimarySources" -> {"P08 Eqs. 27-28", "P09 Eqs. 2.1-2.3"},
   "BLVelocityResiduals" -> <|
    "ut" -> required[[1]], "uphi" -> required[[2]],
    "urSquared" -> urSquaredRes, "urInward" -> inwardRes,
    "utheta" -> 0|>,
   "DefiningEquationResiduals" -> <|
    "E=-u_t" -> energyRes, "Lz=u_phi" -> lzRes,
    "g(u,u)=-1" -> normRes|>,
   "KIdentityResidual" -> kRes,
   "HyperboloidalVelocityResiduals" -> <|
    "uT" -> uTRes, "uR" -> uRRes, "uPhi" -> uPhiRes, "uy" -> 0|>,
   "HorizonRegularChecks" -> <|
    "DeltaAtRPlusResidual" -> deltaPlusRes,
    "RgeoAtRPlusMinusPPlusSquared" -> rGeoPlusRes,
    "ComputedPPlusDefinitionResidual" -> pPlusDefinitionRes,
    "FutureSqrtBranchResidual" -> sqrtBranchRes,
    "KPlusResidual" -> kPlusRes|>,
   "Assumptions" -> {
    "M>0, |a|<M, L>0, r>rPlus, Delta>0",
    "Q=0, theta=Pi/2, Rgeo>=0, inward sign",
    "P+Sqrt[Rgeo]!=0 for K regular form",
    "future BL/hyperboloidal branches require ut>0,uT>0",
    "future horizon branch requires PPlus>0"},
   "FutureDirectionDerivedForAllParameters" -> False,
   "TurningPointDifferentiabilityVerified" -> False,
   "AllRequiredResidualsZero" -> And @@ (# === 0 & /@ required)|>
]
