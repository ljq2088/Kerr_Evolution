(* A4a source blocks: complete point-particle input, compared in jet basis. *)
Module[
 {block, r, del, rho, rhob, ass, extract, l0, j0, lTheta, jTheta,
  base, harmonic, w0, wt, nn0, mn0, mm0, nnt, mnt, mmt,
  generated, target, rawDifference, jetRules, atPoint, jetDifference,
  expandedDifference, numerator, jets, algebraVars, coeffValues,
  coeffResiduals, operandResiduals},
 block = If[ValueQ[a4SourceBlock], a4SourceBlock, 1];
 r = len^2/RR;
 del = r^2 - 2 mass r + spin^2;
 rho = 1/(r - I spin Cos[theta]);
 rhob = 1/(r + I spin Cos[theta]);
 ass = Element[{mass, spin, len, r0, theta0}, Reals] &&
   mass > 0 && len > 0 && r0 > 0 && 0 < theta0 < Pi &&
   (len^2/r0)^2 - 2 mass len^2/r0 + spin^2 != 0 &&
   Element[mmode, Integers];
 extract[z_] := Expand[z/Exp[I mmode PP]];
 l0[z_, s_] := D[z, theta] - I Csc[theta] D[z, PP] -
   I spin Sin[theta] D[z, TT] + s Cot[theta] z;
 j0[z_] := -(2 + 4 mass RR/len^2) D[z, TT] -
   RR^2/len^2 D[z, RR];
 lTheta[z_, s_] := D[z, theta] + mmode Csc[theta] z -
   I spin Sin[theta] D[z, TT] + s Cot[theta] z;
 jTheta[z_] := -(2 + 4 mass RR/len^2) D[z, TT] -
   RR^2/len^2 D[z, RR];
 base = mu len^2/(sigP[TT] rp[TT]^2 Sin[thetaP] uTp[TT])
   deltaR[TT, RR] deltaTheta[theta];
 harmonic = Exp[I mmode (PP - phiP[TT])]/(2 Pi);
 w0 = base harmonic;
 wt = base Exp[-I mmode phiP[TT]]/(2 Pi);
 nn0 = w0 np[TT]^2; mn0 = w0 mbp[TT] np[TT]; mm0 = w0 mbp[TT]^2;
 nnt = wt np[TT]^2; mnt = wt mbp[TT] np[TT]; mmt = wt mbp[TT]^2;
 generated = Switch[block,
   1, extract[-rho^8 rhob l0[rho^-4 l0[rho^-2 rhob^-1 nn0, 0], -1]],
   2, extract[-(1/Sqrt[2]) rho^8 rhob del^2 l0[
      rho^-4 rhob^2 j0[rho^-2 rhob^-2 del^-1 mn0], -1]],
   3, extract[-(1/2) rho^8 rhob del^2 j0[
      rho^-4 j0[rho^-2 rhob mm0]]],
   4, extract[-(1/Sqrt[2]) rho^8 rhob del^2 j0[
      rho^-4 rhob^2 del^-1 l0[rho^-2 rhob^-2 mn0, -1]]]];
 target = Switch[block,
   1, -rho^8 rhob lTheta[rho^-4 lTheta[rho^-2 rhob^-1 nnt, 0], -1],
   2, -(1/Sqrt[2]) rho^8 rhob del^2 lTheta[
      rho^-4 rhob^2 jTheta[rho^-2 rhob^-2 del^-1 mnt], -1],
   3, -(1/2) rho^8 rhob del^2 jTheta[
      rho^-4 jTheta[rho^-2 rhob mmt]],
   4, -(1/Sqrt[2]) rho^8 rhob del^2 jTheta[
      rho^-4 rhob^2 del^-1 lTheta[rho^-2 rhob^-2 mnt, -1]]];
 rawDifference = generated - target;
 jetRules = {
   Derivative[i__][deltaR][x__] :> jetDeltaR[i],
   Derivative[i__][deltaTheta][x__] :> jetDeltaTheta[i],
   Derivative[i__][sigP][x__] :> jetSigP[i],
   Derivative[i__][rp][x__] :> jetRp[i],
   Derivative[i__][uTp][x__] :> jetUTp[i],
   Derivative[i__][phiP][x__] :> jetPhiP[i],
   Derivative[i__][np][x__] :> jetNp[i],
   Derivative[i__][mbp][x__] :> jetMbp[i],
   deltaR[x__] :> jetDeltaR[0], deltaTheta[x__] :> jetDeltaTheta[0],
   sigP[x__] :> jetSigP[0], rp[x__] :> jetRp[0],
   uTp[x__] :> jetUTp[0], phiP[x__] :> jetPhiP[0],
   np[x__] :> jetNp[0], mbp[x__] :> jetMbp[0]};
 atPoint = {TT -> 0, RR -> r0, theta -> theta0, PP -> 0};
 jetDifference = rawDifference /. jetRules /. atPoint;
 expandedDifference = Expand[jetDifference];
 numerator = Numerator[Together[expandedDifference]];
 jets = DeleteDuplicates@Cases[numerator,
   (jetDeltaR | jetDeltaTheta | jetSigP | jetRp | jetUTp |
     jetPhiP | jetNp | jetMbp)[__], Infinity];
 algebraVars = Array[jv, Length[jets]];
 coeffValues = Values@CoefficientRules[
   Expand[numerator /. Thread[jets -> algebraVars]], algebraVars];
 coeffResiduals = FullSimplify[#, ass] & /@ coeffValues;
 operandResiduals = FullSimplify[#, ass] & /@ {
   extract[nn0] - nnt, extract[mn0] - mnt, extract[mm0] - mmt};
 <|"Stage" -> "A4a source jet block pilot",
   "Block" -> block,
   "CompletePointParticleInputDifferentiated" -> True,
   "ArbitraryProjectionTemplateUsed" -> False,
   "OperandModeResiduals" -> AssociationThread[
    {"W*N^2", "W*Mbar*N", "W*Mbar^2"}, operandResiduals],
   "RawDifferenceLeafCount" -> LeafCount[rawDifference],
   "JetBasisSize" -> Length[jets],
   "CoefficientResiduals" -> coeffResiduals,
   "AllCoefficientResidualsZero" -> And @@ (# === 0 & /@ coeffResiduals),
   "ExpandedDifferenceExactlyZero" -> (expandedDifference === 0),
   "Assumptions" -> {
    "M>0,L>0,R0>0,0<theta0<Pi,Delta(R0)!=0,integer m",
    "worldline and delta derivative jets are independent algebraic basis elements"}|>
]
