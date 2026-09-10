(* A6 targeted check: D.4 coefficients and continuous P/Q equivalence. *)
Module[
 {s, d, angular, op, expr, jetRules, jetExpr, extracted, target,
  coeffResiduals, a, b, c, dd, er, f, secondOrder, pDefTime,
  pSystemTime, recoveryResidual, psiTimeFromP, pDefinitionResidual,
  qEvolutionResidual, sourceChecks, ass},
 s = -2; d = 1 - yy^2;
 angular[z_] := d D[z, {yy, 2}] - 2 yy D[z, yy] -
   (mmode - s yy)^2/d z + s z;
 op[z_] :=
   (8 mass (2 mass - spin^2 RR/len^2) (1 + 2 mass RR/len^2) -
      spin^2 d) D[z, {TT, 2}]
   - 2 (len^2 - (8 mass^2 - spin^2) RR^2/len^2 +
      4 spin^2 mass RR^3/len^4) D[z, TT, RR]
   - (len^2 - 2 mass RR + spin^2 RR^2/len^2) RR^2/len^2
      D[z, {RR, 2}] - angular[z]
   + 2 I spin mmode (1 + 4 mass RR/len^2) D[z, TT]
   + 2 I spin mmode RR^2/len^2 D[z, RR]
   + 2 (2 mass (-s + 2 (s + 2) mass RR/len^2 -
      3 spin^2 RR^2/len^4) - spin^2 RR/len^2 -
      I s spin yy) D[z, TT]
   + 2 RR (-(1 + s) + (s + 3) mass RR/len^2 -
      2 spin^2 RR^2/len^4) D[z, RR]
   + 2 I spin mmode RR/len^2 z
   + 2 ((1 + s) mass RR/len^2 - spin^2 RR^2/len^4) z;
 expr = Expand[op[psi[TT, RR, yy]] + angular[psi[TT, RR, yy]]];
 jetRules = {psi[TT, RR, yy] -> j0,
   Derivative[1, 0, 0][psi][TT, RR, yy] -> jT,
   Derivative[0, 1, 0][psi][TT, RR, yy] -> jR,
   Derivative[2, 0, 0][psi][TT, RR, yy] -> jTT,
   Derivative[1, 1, 0][psi][TT, RR, yy] -> jTR,
   Derivative[0, 2, 0][psi][TT, RR, yy] -> jRR};
 jetExpr = expr /. jetRules;
 extracted = AssociationThread[{"A", "B", "C", "D", "ER", "F"},
   Coefficient[jetExpr, #] & /@ {jTT, jTR, jRR, jT, jR, j0}];
 target = <|
   "A" -> 8 mass (2 mass - spin^2 RR/len^2)
      (1 + 2 mass RR/len^2) - spin^2 d,
   "B" -> -2 (len^2 - (8 mass^2 - spin^2) RR^2/len^2 +
      4 spin^2 mass RR^3/len^4),
   "C" -> -(len^2 - 2 mass RR + spin^2 RR^2/len^2) RR^2/len^2,
   "D" -> 2 I spin mmode (1 + 4 mass RR/len^2) +
      2 (2 mass (2 - 3 spin^2 RR^2/len^4) -
       spin^2 RR/len^2 + 2 I spin yy),
   "ER" -> 2 I spin mmode RR^2/len^2 +
      2 RR (1 + mass RR/len^2 - 2 spin^2 RR^2/len^4),
   "F" -> 2 I spin mmode RR/len^2 - 2 mass RR/len^2 -
      2 spin^2 RR^2/len^4|>;
 ass = Element[{mass, spin, len, RR, yy}, Reals] && len > 0 &&
   -1 < yy < 1;
 coeffResiduals = AssociationMap[
   FullSimplify[extracted[#] - target[#], ass] &, Keys[target]];
 a = target["A"]; b = target["B"]; c = target["C"];
 dd = target["D"]; er = target["ER"]; f = target["F"];
 secondOrder = a psiTT + b psiTR + c psiRR - angPsi +
   dd psiT + er psiR + f psi0 - src;
 pDefTime = a psiTT + b psiTR + dd psiT;
 pSystemTime = src - c psiRR + angPsi - er psiR - f psi0;
 recoveryResidual = Expand[pDefTime - pSystemTime - secondOrder];
 psiTimeFromP = (p0 - b q0 - dd psi0)/a;
 pDefinitionResidual = FullSimplify[
   a psiTimeFromP + b q0 + dd psi0 - p0, ass && a != 0];
 qEvolutionResidual = Expand[qT - dRpsiT /. {qT -> psiTR,
    dRpsiT -> psiTR}];
 sourceChecks = <|
   "SourceCoefficientInPT" -> Coefficient[pSystemTime, src],
   "SourceCoefficientInPsiT" -> Coefficient[psiTimeFromP, src],
   "SourceCoefficientInQT" -> Coefficient[qEvolutionResidual, src]|>;
 <|"Stage" -> "A6 first-order identity",
   "ExtractedCoefficientResiduals" -> coeffResiduals,
   "SecondOrderRecoveryResidual" -> recoveryResidual,
   "PDefinitionResidual" -> pDefinitionResidual,
   "QEvolutionIdentityResidual" -> qEvolutionResidual,
   "SourcePlacement" -> sourceChecks,
   "FirstOrderIdentityChecked" -> And @@ (# === 0 & /@
     Join[Values[coeffResiduals], {recoveryResidual,
       pDefinitionResidual, qEvolutionResidual,
       sourceChecks["SourceCoefficientInPT"] - 1,
       sourceChecks["SourceCoefficientInPsiT"],
       sourceChecks["SourceCoefficientInQT"]}]),
   "Assumptions" -> {"L>0", "-1<y<1", "A!=0 for solving psi_T"}|>
]
