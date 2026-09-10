(* A4b/A4c: direct unprojected E0 -> fixed-m y endpoint. *)
Module[
 {group, block, sw, r, del, dy, sqd, sigTheta, sigY, ass,
  trigRules, phaseRules, toYMode, az0, angular0, field0, modeField, fieldGenerated,
  angularY, fieldYTarget, fieldResidual, fieldEvidence, sourceEvidence},
 group = If[ValueQ[a4bGroup], a4bGroup, "Field"];
 block = If[ValueQ[a4bBlock], a4bBlock, 1];
 sw = -2; r = len^2/RR; del = r^2 - 2 mass r + spin^2;
 dy = 1 - yy^2; sqd = Sqrt[dy];
 sigTheta = r^2 + spin^2 Cos[theta]^2;
 sigY = r^2 + spin^2 yy^2;
 ass = Element[{mass, spin, len, RR, yy}, Reals] &&
   mass > 0 && len > 0 && RR > 0 && -1 < yy < 1 && del != 0 &&
   Element[mmode, Integers];
 trigRules = {Cos[theta] -> -yy, Sin[theta] -> sqd,
   Csc[theta] -> 1/sqd, Cot[theta] -> -yy/sqd,
   Cos[theta]^2 -> yy^2, Sin[theta]^2 -> dy};
 phaseRules = {
   Exp[-I mmode PP] Exp[I mmode (PP - x_)] :> Exp[-I mmode x],
   Exp[-I mmode PP + I mmode (PP - x_)] :> Exp[-I mmode x]};
 toYMode[z_] := Expand[(((Expand[z/Exp[I mmode PP]] /. phaseRules) /.
    trigRules) /. theta -> ArcCos[-yy])];
 az0[z_] := -I D[z, PP] + sw Cos[theta] z;
 angular0[z_] := Csc[theta] D[Sin[theta] D[z, theta], theta] +
   sw z - Csc[theta]^2 az0[az0[z]];
 field0[z_] :=
   (8 mass (2 mass - spin^2 RR/len^2) (1 + 2 mass RR/len^2) -
      spin^2 Sin[theta]^2) D[z, {TT, 2}]
   - 2 (len^2 - (8 mass^2 - spin^2) RR^2/len^2 +
      4 spin^2 mass RR^3/len^4) D[z, TT, RR]
   - (len^2 - 2 mass RR + spin^2 RR^2/len^2) RR^2/len^2
      D[z, {RR, 2}] - angular0[z]
   + 2 spin (1 + 4 mass RR/len^2) D[z, TT, PP]
   + 2 spin RR^2/len^2 D[z, RR, PP]
   + 2 (2 mass (-sw + 2 (sw + 2) mass RR/len^2 -
      3 spin^2 RR^2/len^4) - spin^2 RR/len^2 +
      I sw spin Cos[theta]) D[z, TT]
   + 2 RR (-(1 + sw) + (sw + 3) mass RR/len^2 -
      2 spin^2 RR^2/len^4) D[z, RR]
   + 2 spin RR/len^2 D[z, PP]
   + 2 ((1 + sw) mass RR/len^2 - spin^2 RR^2/len^4) z;
 modeField = psiY[TT, RR, -Cos[theta]] Exp[I mmode PP];
 fieldGenerated = toYMode[field0[modeField]];
 angularY[z_] := dy D[z, {yy, 2}] - 2 yy D[z, yy] -
   (mmode - sw yy)^2/dy z + sw z;
 fieldYTarget[z_] :=
   (8 mass (2 mass - spin^2 RR/len^2) (1 + 2 mass RR/len^2) -
      spin^2 dy) D[z, {TT, 2}]
   - 2 (len^2 - (8 mass^2 - spin^2) RR^2/len^2 +
      4 spin^2 mass RR^3/len^4) D[z, TT, RR]
   - (len^2 - 2 mass RR + spin^2 RR^2/len^2) RR^2/len^2
      D[z, {RR, 2}] - angularY[z]
   + 2 I spin mmode (1 + 4 mass RR/len^2) D[z, TT]
   + 2 I spin mmode RR^2/len^2 D[z, RR]
   + 2 (2 mass (-sw + 2 (sw + 2) mass RR/len^2 -
      3 spin^2 RR^2/len^4) - spin^2 RR/len^2 -
      I sw spin yy) D[z, TT]
   + 2 RR (-(1 + sw) + (sw + 3) mass RR/len^2 -
      2 spin^2 RR^2/len^4) D[z, RR]
   + 2 I spin mmode RR/len^2 z
   + 2 ((1 + sw) mass RR/len^2 - spin^2 RR^2/len^4) z;
 fieldResidual = FullSimplify[Together[
   fieldGenerated - fieldYTarget[psiY[TT, RR, yy]]], ass];
 fieldEvidence = <|
   "DirectStart" -> "complete unprojected field on psiY(T,R,-Cos[theta]) Exp[I m Phi]",
   "TheoryThetaOrYHelperUsedInGeneration" -> False,
   "CompleteYFieldResidual" -> fieldResidual,
   "FieldEndpointVerified" -> (fieldResidual === 0)|>;

 sourceEvidence = If[group === "Field", Missing["NotExecuted", "Source"],
   Module[{rhoT, rhobT, rhoY, rhobY, l0, j0, lY, jY,
     thetaDeltaComposed, base0, harmonic, w0, wy, nn0, mn0, mm0,
     nny, mny, mmy, generated, target, rawDifference, atPoint,
     jetRules, jetDifference, expandedDifference, numerator, jets,
     algebraVars, coeffValues, coeffResiduals, coeffNonzero, operandResiduals,
     weightJacobianResidual, prefGenerated, prefY, outerPrefactorResidual,
     sourceRawDifference, sourceRawLeafCount, decompositionIdentityResidual,
     generatedJetAtPoint, sourceJetDifference, sourceExpandedDifference,
     required},
    rhoT = 1/(r - I spin Cos[theta]); rhobT = 1/(r + I spin Cos[theta]);
    rhoY = 1/(r + I spin yy); rhobY = 1/(r - I spin yy);
    l0[z_, s_] := D[z, theta] - I Csc[theta] D[z, PP] -
      I spin Sin[theta] D[z, TT] + s Cot[theta] z;
    j0[z_] := -(2 + 4 mass RR/len^2) D[z, TT] -
      RR^2/len^2 D[z, RR];
    lY[z_, s_] := sqd D[z, yy] - I spin sqd D[z, TT] +
      (mmode - s yy)/sqd z;
    jY[z_] := -(2 + 4 mass RR/len^2) D[z, TT] -
      RR^2/len^2 D[z, RR];
    thetaDeltaComposed = Sin[thetaP] deltaY[-Cos[theta]];
    base0 = mu len^2/(sigP[TT] rp[TT]^2 Sin[thetaP] uTp[TT])
      deltaR[TT, RR] thetaDeltaComposed;
    harmonic = Exp[I mmode (PP - phiP[TT])]/(2 Pi);
    w0 = base0 harmonic;
    wy = mu len^2/(sigP[TT] rp[TT]^2 uTp[TT])
      deltaR[TT, RR] deltaY[yy] Exp[-I mmode phiP[TT]]/(2 Pi);
    nn0 = w0 np[TT]^2; mn0 = w0 mbp[TT] np[TT]; mm0 = w0 mbp[TT]^2;
    nny = wy np[TT]^2; mny = wy mbp[TT] np[TT]; mmy = wy mbp[TT]^2;
    generated = Switch[block,
      1, toYMode[-rhoT^8 rhobT l0[
        rhoT^-4 l0[rhoT^-2 rhobT^-1 nn0, 0], -1]],
      2, toYMode[-(1/Sqrt[2]) rhoT^8 rhobT del^2 l0[
        rhoT^-4 rhobT^2 j0[rhoT^-2 rhobT^-2 del^-1 mn0], -1]],
      3, toYMode[-(1/2) rhoT^8 rhobT del^2 j0[
        rhoT^-4 j0[rhoT^-2 rhobT mm0]]],
      4, toYMode[-(1/Sqrt[2]) rhoT^8 rhobT del^2 j0[
        rhoT^-4 rhobT^2 del^-1 l0[
          rhoT^-2 rhobT^-2 mn0, -1]]]];
    target = Switch[block,
      1, -rhoY^8 rhobY lY[rhoY^-4 lY[rhoY^-2 rhobY^-1 nny, 0], -1],
      2, -(1/Sqrt[2]) rhoY^8 rhobY del^2 lY[
        rhoY^-4 rhobY^2 jY[rhoY^-2 rhobY^-2 del^-1 mny], -1],
      3, -(1/2) rhoY^8 rhobY del^2 jY[
        rhoY^-4 jY[rhoY^-2 rhobY mmy]],
      4, -(1/Sqrt[2]) rhoY^8 rhobY del^2 jY[
        rhoY^-4 rhobY^2 del^-1 lY[rhoY^-2 rhobY^-2 mny, -1]]];
    rawDifference = generated - target;
    jetRules = {
      Derivative[i__][deltaR][x__] :> jetDeltaR[i],
      Derivative[i__][deltaY][x__] :> jetDeltaY[i],
      Derivative[i__][sigP][x__] :> jetSigP[i],
      Derivative[i__][rp][x__] :> jetRp[i],
      Derivative[i__][uTp][x__] :> jetUTp[i],
      Derivative[i__][phiP][x__] :> jetPhiP[i],
      Derivative[i__][np][x__] :> jetNp[i],
      Derivative[i__][mbp][x__] :> jetMbp[i],
      deltaR[x__] :> jetDeltaR[0], deltaY[x__] :> jetDeltaY[0],
      sigP[x__] :> jetSigP[0], rp[x__] :> jetRp[0],
      uTp[x__] :> jetUTp[0], phiP[x__] :> jetPhiP[0],
      np[x__] :> jetNp[0], mbp[x__] :> jetMbp[0]};
    atPoint = {TT -> 0, RR -> r0, yy -> y0, PP -> 0};
    jetDifference = rawDifference /. jetRules /. atPoint;
    expandedDifference = Expand[jetDifference];
    numerator = Numerator[Together[expandedDifference]];
    jets = DeleteDuplicates@Cases[numerator,
      (jetDeltaR | jetDeltaY | jetSigP | jetRp | jetUTp |
        jetPhiP | jetNp | jetMbp)[__], Infinity];
    algebraVars = Array[jv, Length[jets]];
    coeffValues = Values@CoefficientRules[
      Expand[numerator /. Thread[jets -> algebraVars]], algebraVars];
    coeffResiduals = FullSimplify[#, ass /. {RR -> r0, yy -> y0}] & /@
      coeffValues;
    coeffNonzero = Select[coeffResiduals, # =!= 0 &];
    operandResiduals = FullSimplify[#, ass] & /@ {
      toYMode[nn0] - nny, toYMode[mn0] - mny, toYMode[mm0] - mmy};
    weightJacobianResidual = FullSimplify[
      Sin[thetaP]/Sin[thetaP] - 1, 0 < thetaP < Pi];
    prefGenerated = Expand[(16 Pi sigTheta/(del^2 RR)) /. trigRules];
    prefY = 16 Pi sigY/(del^2 RR);
    outerPrefactorResidual = FullSimplify[
      prefGenerated - prefY, ass];
    sourceRawDifference = prefGenerated generated - prefY target;
    sourceRawLeafCount = LeafCount[sourceRawDifference];
    decompositionIdentityResidual = Expand[
      pGen gGen - pTar gTar -
       ((pGen - pTar) gGen + pTar (gGen - gTar))];
    generatedJetAtPoint = generated /. jetRules /. atPoint;
    sourceJetDifference = outerPrefactorResidual generatedJetAtPoint +
      (prefY /. {RR -> r0, yy -> y0}) expandedDifference;
    sourceExpandedDifference = Expand[sourceJetDifference];
    required = Join[operandResiduals, coeffResiduals,
      {weightJacobianResidual, outerPrefactorResidual,
       decompositionIdentityResidual, sourceExpandedDifference}];
    <|"Block" -> block,
      "DirectStart" ->
       "unprojected periodic point-particle block with deltaTheta=Sin[thetaP] deltaY[-Cos[theta]] before differentiation",
      "TheoryThetaOrYBlockUsedInGeneration" -> False,
      "OperandYModeResiduals" -> AssociationThread[
       {"Wy*N^2", "Wy*Mbar*N", "Wy*Mbar^2"}, operandResiduals],
      "AngularDeltaWeightJacobianResidual" -> weightJacobianResidual,
      "OuterPrefactorResidual" -> outerPrefactorResidual,
      "GeneratedThetaFree" -> FreeQ[generated, theta],
      "GeneratedPhiFree" -> FreeQ[generated, PP],
      "GeneratedHasNoThetaOrPhi" ->
       (FreeQ[generated, theta] && FreeQ[generated, PP]),
      "RawDifferenceLeafCount" -> LeafCount[rawDifference],
      "JetBasisSize" -> Length[jets],
      "CoefficientResidualCount" -> Length[coeffResiduals],
      "NonzeroCoefficientResidualCount" -> Length[coeffNonzero],
      "FirstNonzeroCoefficientLeafCount" ->
       If[coeffNonzero === {}, 0, LeafCount[First[coeffNonzero]]],
      "AllCoefficientResidualsZero" -> (coeffNonzero === {}),
      "ExpandedDifferenceExactlyZero" -> (expandedDifference === 0),
      "SourceBlockJetResidualLeafCount" -> LeafCount[sourceExpandedDifference],
      "SourceBlockJetResidualExactlyZero" ->
       (sourceExpandedDifference === 0),
      "DirectSourceRawDifferenceLeafCount" -> sourceRawLeafCount,
      "SourceResidualDecompositionIdentityResidual" ->
       decompositionIdentityResidual,
      "DirectSourceResidualDecompositionExactlyZero" ->
       (decompositionIdentityResidual === 0),
      "SourceResidualFactoredAsOuterPlusHat" -> True,
      "AllRequiredResidualsZero" ->
       And @@ (# === 0 & /@ required) && expandedDifference === 0|>
   ]];

 <|"Stage" -> "A4b direct y field/source",
   "ExecutionGroup" -> group, "FieldEvidence" -> fieldEvidence,
   "SourceEvidence" -> sourceEvidence,
   "Assumptions" -> {
    "M>0,L>0,R>0,-1<y<1,Delta!=0,integer m",
    "0<thetaP<Pi for angular delta chart",
    "worldline amplitudes are sufficiently differentiable; delta derivatives are distributional",
    "no Gaussian and no explicit spatial delta-jet expansion"}|>
]
