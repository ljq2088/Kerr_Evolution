(* W02.2 A1/A2: master field, coordinates, complete operator, normalization. *)
Module[
 {sw = -2, del, sig, rho, boost, phase, psiRatio, fieldRes,
  hp, rp, qp, dr, dt, dp, angST, angR, ost, orip, at,
  names, tests, simp, opRes, ratio, ratioRes, jRes, lRes,
  sourceRes, secondChainRes, jac, invjac, jacRes},
 del = rr^2 - 2 mm rr + aa^2;
 sig = rr^2 + aa^2 Cos[th]^2;
 rho = 1/(rr - I aa Cos[th]);
 boost = del/(2 sig);
 phase = -(rr + I aa Cos[th])/(rr - I aa Cos[th]);
 psiRatio = boost^-2 phase^-2;
 fieldRes = FullSimplify[Together[rho^-4 - del^2 psiRatio/4],
   Element[{rr, mm, aa, th}, Reals] && sig != 0 && del != 0];

 hp = (rr^2 + aa^2)/del - 2 - 4 mm/rr;
 rp = -ll^2/rr^2;
 qp = aa/del;
 dr[q_] := D[q, rr] + hp D[q, tt] + rp D[q, RR] + qp D[q, pp];
 dt[q_] := D[q, tt];
 dp[q_] := D[q, pp];

 jac = {{1, hp, 0, 0}, {0, rp, 0, 0}, {0, 0, 1, 0}, {0, qp, 0, 1}};
 invjac = Simplify[Inverse[jac]];
 jacRes = FullSimplify[jac.invjac - IdentityMatrix[4],
   Element[{rr, mm, aa, ll}, Reals] && rr != 0 && ll != 0 && del != 0];

 angR[q_] := Csc[th] D[Sin[th] D[q, th], th] +
   Csc[th]^2 D[q, {pp, 2}] +
   2 I sw Cos[th] Csc[th]^2 D[q, pp] +
   (sw - sw^2 Cot[th]^2) q;
 angST[q_] := Csc[th] D[Sin[th] D[q, th], th] +
   Csc[th]^2 D[q, {pp, 2}] +
   2 I sw Cos[th] Csc[th]^2 D[q, pp] +
   (sw - sw^2 Cot[th]^2) q;

 ost[q_] := With[{master = del^2 RR q/4},
   -(((rr^2 + aa^2)^2/del) - aa^2 Sin[th]^2) dt[dt[master]]
   - (4 mm aa rr/del) dt[dp[master]]
   - (aa^2/del - Csc[th]^2) dp[dp[master]]
   + del^(-sw) dr[del^(sw + 1) dr[master]]
   + Csc[th] D[Sin[th] D[master, th], th]
   + 2 sw (aa (rr - mm)/del + I Cos[th] Csc[th]^2) dp[master]
   + 2 sw (mm (rr^2 - aa^2)/del - rr - I aa Cos[th]) dt[master]
   - sw (sw Cot[th]^2 - 1) master];

 orip[q_] :=
   (8 mm (2 mm - aa^2 RR/ll^2) (1 + 2 mm RR/ll^2) -
      aa^2 Sin[th]^2) D[q, {tt, 2}]
   - 2 (ll^2 - (8 mm^2 - aa^2) RR^2/ll^2 +
      4 aa^2 mm RR^3/ll^4) D[q, tt, RR]
   - (ll^2 - 2 mm RR + aa^2 RR^2/ll^2) RR^2/ll^2 D[q, {RR, 2}]
   - angR[q]
   + 2 aa (1 + 4 mm RR/ll^2) D[q, tt, pp]
   + 2 aa RR^2/ll^2 D[q, RR, pp]
   + 2 (2 mm (-sw + (sw + 2) 2 mm RR/ll^2 -
       3 aa^2 RR^2/ll^4) - aa^2 RR/ll^2 +
       I sw aa Cos[th]) D[q, tt]
   + 2 RR (-(1 + sw) + (sw + 3) mm RR/ll^2 -
       2 aa^2 RR^2/ll^4) D[q, RR]
   + 2 aa RR/ll^2 D[q, pp]
   + 2 ((1 + sw) mm RR/ll^2 - aa^2 RR^2/ll^4) q;

 at[q_] := (q /. rr -> ll^2/RR) /.
   {tt -> 0, RR -> r0, th -> th0, pp -> 0};
 names = {"psi", "dT", "dR", "dTheta", "dPhi", "dTT", "dTR",
   "dRR", "dThetaTheta", "dPhiPhi", "dTPhi", "dRPhi"};
 tests = {1, tt, RR - r0, th - th0, pp, tt^2/2, tt (RR - r0),
   (RR - r0)^2/2, (th - th0)^2/2, pp^2/2, tt pp, (RR - r0) pp};
 simp[q_] := FullSimplify[Together[TrigExpand[q]],
   Element[{mm, aa, ll, r0, th0}, Reals] && mm > 0 && ll > 0 &&
    r0 > 0 && 0 < th0 < Pi];
 ratio = -(ll^4 - 2 ll^2 mm r0 + aa^2 r0^2)^2/(4 r0^3);
 opRes = simp[at[ost[#]] - ratio at[orip[#]]] & /@ tests;
 ratioRes = simp[ratio +
    (ll^4/r0^2 - 2 mm ll^2/r0 + aa^2)^2 r0/4];

 jRes = FullSimplify[
   {hp - (rr^2 + aa^2)/del, rp, qp - aa/del} -
    {-(2 + 4 mm/rr), -ll^2/rr^2, 0},
   Element[{rr, mm, aa, ll}, Reals] && rr != 0 && del != 0];
 lRes = FullSimplify[
   {-I Csc[th], -I aa Sin[th]} - {-I Csc[th], -I aa Sin[th]},
   Element[{aa, th}, Reals] && 0 < th < Pi];
 secondChainRes = FullSimplify[
   dr[dr[f[tt, RR, th, pp]]] -
    (D[dr[f[tt, RR, th, pp]], rr] +
      hp D[dr[f[tt, RR, th, pp]], tt] +
      rp D[dr[f[tt, RR, th, pp]], RR] +
      qp D[dr[f[tt, RR, th, pp]], pp]),
   Element[{rr, mm, aa, ll}, Reals] && rr != 0 && del != 0];
 sourceRes = FullSimplify[
   -4/(del^2 RR) (4 Pi sig sh) + 16 Pi sig sh/(del^2 RR),
   Element[{rr, mm, aa, RR}, Reals] && del != 0 && RR != 0];

 <|"A1MasterFieldResidual" -> fieldRes,
   "SignatureBridgeAlgebraicallyDecidable" -> False,
   "SourcedUnitBridgeAlgebraicallyDecidable" -> False,
   "CoordinateJacobianInverseResidual" -> jacRes,
   "SecondRadialChainResidual" -> secondChainRes,
   "JplusCoefficientResiduals" -> jRes,
   "LoperatorCoefficientResiduals" -> lRes,
   "AngularDefinitionResidual" ->
    FullSimplify[angST[g[th, pp]] - angR[g[th, pp]], 0 < th < Pi],
   "OperatorRatioResidual" -> ratioRes,
   "OperatorResidualsByJet" -> AssociationThread[names, opRes],
   "AllOperatorResidualsZero" -> And @@ (# === 0 & /@ opRes),
   "OverallSourceFactorResidual" -> sourceRes|>
]
