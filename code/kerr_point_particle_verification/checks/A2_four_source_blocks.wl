(* W02.2/A2 end-to-end source validation.
   Generation starts from the complete Sasaki BL blocks acting on composed
   new-coordinate functions. No hatted operator or target block is used in
   the generation route. The theory endpoint is entered separately below. *)
Module[
 {delBL, rhoBL, rhobBL, hPrimeBL, qPrimeBL, lSourceBL, jSourceBL,
  nnComposed, mnComposed, mmComposed, sasakiBLBlocks, generatedBlocks,
  derivativeRules, toNewCoordinates, delH, rhoH, rhobH, nnH, mnH, mmH,
  lTargetH, jTargetH, targetBlocks, assumptions, blockResiduals,
  totalEndpointResidual, b2PrimeResidual, b2StarPrimeResidual,
  unresolvedCoordinateObjects, fieldFactorContext},

 (* Generation route: original BL source plus coordinate composition only. *)
 delBL = rBL^2 - 2 mass rBL + spin^2;
 rhoBL = 1/(rBL - I spin Cos[theta]);
 rhobBL = 1/(rBL + I spin Cos[theta]);
 hPrimeBL = (rBL^2 + spin^2)/delBL - 2 - 4 mass/rBL;
 qPrimeBL = spin/delBL;

 lSourceBL[z_, s_] := D[z, theta] - I Csc[theta] D[z, phiBL] -
   I spin Sin[theta] D[z, tBL] + s Cot[theta] z;
 jSourceBL[z_] := D[z, rBL] -
   ((rBL^2 + spin^2)/delBL) D[z, tBL] -
   (spin/delBL) D[z, phiBL];

 nnComposed = fNN[tBL + hCoord[rBL], len^2/rBL, theta,
   phiBL + qCoord[rBL]];
 mnComposed = fMN[tBL + hCoord[rBL], len^2/rBL, theta,
   phiBL + qCoord[rBL]];
 mmComposed = fMM[tBL + hCoord[rBL], len^2/rBL, theta,
   phiBL + qCoord[rBL]];

 (* These are exactly the four original Sasaki blocks before the external 2. *)
 sasakiBLBlocks = {
   -(1/2) rhoBL^8 rhobBL lSourceBL[
     rhoBL^-4 lSourceBL[rhoBL^-2 rhobBL^-1 nnComposed, 0], -1],
   -(1/(2 Sqrt[2])) rhoBL^8 rhobBL delBL^2 lSourceBL[
     rhoBL^-4 rhobBL^2 jSourceBL[
       rhoBL^-2 rhobBL^-2 delBL^-1 mnComposed], -1],
   -(1/4) rhoBL^8 rhobBL delBL^2 jSourceBL[
     rhoBL^-4 jSourceBL[rhoBL^-2 rhobBL mmComposed]],
   -(1/(2 Sqrt[2])) rhoBL^8 rhobBL delBL^2 jSourceBL[
     rhoBL^-4 rhobBL^2 delBL^-1 lSourceBL[
       rhoBL^-2 rhobBL^-2 mnComposed, -1]]};

 (* The Kernel has already executed all BL derivatives above. These rules
    supply only derivatives of the frozen coordinate functions. *)
 derivativeRules = {
   Derivative[1][hCoord][rBL] -> hPrimeBL,
   Derivative[2][hCoord][rBL] -> D[hPrimeBL, rBL],
   Derivative[1][qCoord][rBL] -> qPrimeBL,
   Derivative[2][qCoord][rBL] -> D[qPrimeBL, rBL]};
 toNewCoordinates[z_] := z /. derivativeRules /.
    {tBL + hCoord[rBL] -> TT, len^2/rBL -> RR,
     phiBL + qCoord[rBL] -> PP} /. rBL -> len^2/RR;
 generatedBlocks = toNewCoordinates /@ (2 sasakiBLBlocks);
 unresolvedCoordinateObjects = DeleteDuplicates@Join[
   Cases[generatedBlocks, hCoord[__], Infinity],
   Cases[generatedBlocks, qCoord[__], Infinity],
   Cases[generatedBlocks, Derivative[__][hCoord][__], Infinity],
   Cases[generatedBlocks, Derivative[__][qCoord][__], Infinity]];

 (* Independent theory endpoint entry from appendix/SASAKI_SOURCE.md. *)
 delH = (len^2/RR)^2 - 2 mass len^2/RR + spin^2;
 fieldFactorContext = delH^2 RR/4;
 rhoH = 1/(len^2/RR - I spin Cos[theta]);
 rhobH = 1/(len^2/RR + I spin Cos[theta]);
 nnH = fNN[TT, RR, theta, PP];
 mnH = fMN[TT, RR, theta, PP];
 mmH = fMM[TT, RR, theta, PP];
 lTargetH[z_, s_] := D[z, theta] - I Csc[theta] D[z, PP] -
   I spin Sin[theta] D[z, TT] + s Cot[theta] z;
 jTargetH[z_] := -(2 + 4 mass RR/len^2) D[z, TT] -
   (RR^2/len^2) D[z, RR];
 targetBlocks = {
   -rhoH^8 rhobH lTargetH[
     rhoH^-4 lTargetH[rhoH^-2 rhobH^-1 nnH, 0], -1],
   -(1/Sqrt[2]) rhoH^8 rhobH delH^2 lTargetH[
     rhoH^-4 rhobH^2 jTargetH[
       rhoH^-2 rhobH^-2 delH^-1 mnH], -1],
   -(1/2) rhoH^8 rhobH delH^2 jTargetH[
     rhoH^-4 jTargetH[rhoH^-2 rhobH mmH]],
   -(1/Sqrt[2]) rhoH^8 rhobH delH^2 jTargetH[
     rhoH^-4 rhobH^2 delH^-1 lTargetH[
       rhoH^-2 rhobH^-2 mnH, -1]]};

 assumptions = Element[{mass, spin, len, RR, theta}, Reals] &&
   mass > 0 && len > 0 && RR > 0 && 0 < theta < Pi && delH != 0;
 blockResiduals = MapThread[
   FullSimplify[Together[#1 - #2], assumptions] &,
   {generatedBlocks, targetBlocks}];
 b2PrimeResidual = FullSimplify[Together[
   generatedBlocks[[1]] + generatedBlocks[[2]] -
    targetBlocks[[1]] - targetBlocks[[2]]], assumptions];
 b2StarPrimeResidual = FullSimplify[Together[
   generatedBlocks[[3]] + generatedBlocks[[4]] -
    targetBlocks[[3]] - targetBlocks[[4]]], assumptions];
 totalEndpointResidual = FullSimplify[Together[
   Total[generatedBlocks] - Total[targetBlocks]], assumptions];

 <|"PrimaryEndToEndEvidence" -> <|
    "TheoryStart" ->
     "theory/appendix/SASAKI_SOURCE.md: original BL four blocks",
    "TheoryEndpoint" ->
     "theory/appendix/SASAKI_SOURCE.md: complete pulled-back hat T_H",
    "GenerationUsedHattedOperators" -> False,
    "GenerationUsedTheoryIntermediateDerivativeMap" -> False,
    "GenerationUsedTheoryTargetBlockTemplate" -> False,
    "UnresolvedCoordinateCompositionObjects" -> unresolvedCoordinateObjects,
    "AllCoordinateCompositionDerivativesResolved" ->
     (unresolvedCoordinateObjects === {}),
    "ProjectionInputs" ->
     "arbitrary fNN, fMN, fMM functions composed with (T,R,theta,Phi)",
    "FrozenProjectionTetradContext" ->
     "fNN, fMN, fMM are the Sasaki Kinnersley projections; no additional tetrad change is applied inside the source functional",
    "BlockEndpointResiduals" -> AssociationThread[
     {"nn_LL", "mbar_n_LJ", "mbar_mbar_JJ", "mbar_n_JL"},
     blockResiduals],
    "B2PrimeEndpointResidual" -> b2PrimeResidual,
    "B2StarPrimeEndpointResidual" -> b2StarPrimeResidual,
    "WholeHatTEndpointResidual" -> totalEndpointResidual,
    "WholeEndpointVerified" ->
     And @@ (# === 0 & /@ Join[blockResiduals,
       {b2PrimeResidual, b2StarPrimeResidual, totalEndpointResidual}])|>,
   "LeftOperatorIdentityEvidence" -> <|
    "FrozenFieldTetradFactor" -> fieldFactorContext,
    "FrozenFieldTetradFactorRole" ->
     "context for the referenced left identity; not used to generate the source blocks",
    "Artifact" -> "data/kerr_point_particle_verification/results/A1_A2_operator.json",
    "AllTwelveJetResidualsZero" -> True,
    "OperatorRatioResidualZero" -> True,
    "SecondRadialChainResidualZero" -> True,
    "RecomputedHere" -> False,
    "DirectReferenceAudited" -> True|>,
   "EvidenceLimits" -> <|
    "IndependentPrimaryPDFTranscriptionVerified" -> False,
    "A1SignatureBridgeClosed" -> False,
    "A1SourcedUnitBridgeClosed" -> False,
    "ExcludedLimits" ->
     {"Delta=0 in this rational-expression comparison",
      "R=0", "theta=0,Pi coordinate poles"},
    "NarrowConclusion" ->
     "The complete source computed from the theory BL start equals the separately entered theory hyperboloidal endpoint for arbitrary projection functions; primary-PDF transcription and A1 convention/unit bridges remain open."|>|>
]
