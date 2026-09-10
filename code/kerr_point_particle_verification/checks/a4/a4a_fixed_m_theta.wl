(* A4a complete-field check only. Source blocks use a4a_source_jet_blocks.wl. *)
Module[
 {sw, r, del, ass, modeField, az0, angular0, field0, extractMode,
  fieldGenerated, angularThetaTarget, fieldThetaTarget, fieldResidual},
 sw = -2;
 r = len^2/RR;
 del = r^2 - 2 mass r + spin^2;
 ass = Element[{mass, spin, len, RR, theta}, Reals] &&
   mass > 0 && len > 0 && RR > 0 && 0 < theta < Pi && del != 0 &&
   Element[mmode, Integers];
 modeField = psiM[TT, RR, theta] Exp[I mmode PP];
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
 extractMode[z_] := FullSimplify[z/Exp[I mmode PP], ass];
 fieldGenerated = extractMode[field0[modeField]];

 angularThetaTarget[z_] := Csc[theta] D[Sin[theta] D[z, theta], theta] +
   (sw - (mmode + sw Cos[theta])^2 Csc[theta]^2) z;
 fieldThetaTarget[z_] :=
   (8 mass (2 mass - spin^2 RR/len^2) (1 + 2 mass RR/len^2) -
      spin^2 Sin[theta]^2) D[z, {TT, 2}]
   - 2 (len^2 - (8 mass^2 - spin^2) RR^2/len^2 +
      4 spin^2 mass RR^3/len^4) D[z, TT, RR]
   - (len^2 - 2 mass RR + spin^2 RR^2/len^2) RR^2/len^2
      D[z, {RR, 2}] - angularThetaTarget[z]
   + 2 I spin mmode (1 + 4 mass RR/len^2) D[z, TT]
   + 2 I spin mmode RR^2/len^2 D[z, RR]
   + 2 (2 mass (-sw + 2 (sw + 2) mass RR/len^2 -
      3 spin^2 RR^2/len^4) - spin^2 RR/len^2 +
      I sw spin Cos[theta]) D[z, TT]
   + 2 RR (-(1 + sw) + (sw + 3) mass RR/len^2 -
      2 spin^2 RR^2/len^4) D[z, RR]
   + 2 I spin mmode RR/len^2 z
   + 2 ((1 + sw) mass RR/len^2 - spin^2 RR^2/len^4) z;
 fieldResidual = FullSimplify[Together[
   fieldGenerated - fieldThetaTarget[psiM[TT, RR, theta]]], ass];
 <|"Stage" -> "A4a fixed-m theta complete field",
   "GenerationStart" ->
    "complete unprojected Ripley operator on psiM(T,R,theta) Exp[I m Phi]",
   "TheoryThetaHelperUsedInGeneration" -> False,
   "CompleteFieldResidual" -> fieldResidual,
   "GeneratedExpressionHasNoPhi" -> FreeQ[fieldGenerated, PP],
   "FieldEndpointVerified" -> (fieldResidual === 0),
   "SourceBlockEvidence" -> "checks/a4/a4a_source_jet_blocks.wl",
   "SourceBlocksExecutedSeparately" -> True,
   "ContainsSourceFullSimplifyBranch" -> False|>
]
