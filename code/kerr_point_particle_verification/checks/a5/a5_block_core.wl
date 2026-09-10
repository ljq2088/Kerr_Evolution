(* Shared executable core for A5b--A5e. Wrappers set a5Block. *)
Module[
 {block, hd, hs, cY, ell, jT, jR, rho, rhob, delta, extPref,
  lDirect, jDirect, fExt, hatBlock, generatedRaw, timeRules,
  generated, f0, f1, f2, lJet, jJet, target, u0, u1, u2,
  v0, v1, w0, w1, residual, expandedResidual, timeRemnants,
  undefinedJetRemnants, ordering, ass},
 block = If[ValueQ[a5Block], a5Block, 1];
 hd = len^4 - 2 mass len^2 RR + spin^2 RR^2;
 hs = len^4 + spin^2 RR^2 yy^2;
 cY = Sqrt[1 - yy^2];
 ell[s_] := (mmode - s yy)/cY;
 jT = -(2 + 4 mass RR/len^2);
 jR = -RR^2/len^2;
 rho = RR/(len^2 + I spin RR yy);
 rhob = RR/(len^2 - I spin RR yy);
 delta = hd/RR^2;
 extPref = -16 Pi RR hs/hd^2;
 lDirect[z_, s_] := cY D[z, yy] - I spin cY D[z, TT] + ell[s] z;
 jDirect[z_] := jT D[z, TT] + jR D[z, RR];
 fExt = Switch[block, 1, hd^2 ff[TT, RR, yy],
   2 | 4, hd ff[TT, RR, yy], 3, ff[TT, RR, yy]];
 hatBlock = Switch[block,
   1, -rho^8 rhob lDirect[
      rho^-4 lDirect[rho^-2 rhob^-1 fExt, 0], -1],
   2, -(1/Sqrt[2]) rho^8 rhob delta^2 lDirect[
      rho^-4 rhob^2 jDirect[rho^-2 rhob^-2 delta^-1 fExt], -1],
   3, -(1/2) rho^8 rhob delta^2 jDirect[
      rho^-4 jDirect[rho^-2 rhob fExt]],
   4, -(1/Sqrt[2]) rho^8 rhob delta^2 jDirect[
      rho^-4 rhob^2 delta^-1 lDirect[
       rho^-2 rhob^-2 fExt, -1]]];
 generatedRaw = extPref hatBlock;
 f0 = fJet0[RR, yy]; f1 = fJet1[RR, yy]; f2 = fJet2[RR, yy];
 timeRules = {
   Derivative[kt_, kr_, ky_][ff][TT, RR, yy] :>
    D[Switch[kt, 0, fJet0[RR, yy], 1, fJet1[RR, yy],
       2, fJet2[RR, yy]], {RR, kr}, {yy, ky}],
   ff[TT, RR, yy] :> fJet0[RR, yy]};
 generated = generatedRaw /. timeRules;

 (* Appendix E.8 endpoint, entered independently from generation. *)
 lJet[s_, xk_, xkp1_] := cY D[xk, yy] - I spin cY xkp1 + ell[s] xk;
 jJet[xk_, xkp1_] := jT xkp1 + jR D[xk, RR];
 Switch[block,
  1,
   u0 = rho^-2 rhob^-1 f0; u1 = rho^-2 rhob^-1 f1;
   u2 = rho^-2 rhob^-1 f2;
   v0 = lJet[0, u0, u1]; v1 = lJet[0, u1, u2];
   w0 = rho^-4 v0; w1 = rho^-4 v1;
   target = -16 Pi RR hs (-rho^8 rhob lJet[-1, w0, w1]);
   ordering = {"inner L_0", "outer L_-1"},
  2,
   u0 = rho^-2 rhob^-2 RR^2 f0; u1 = rho^-2 rhob^-2 RR^2 f1;
   u2 = rho^-2 rhob^-2 RR^2 f2;
   v0 = jJet[u0, u1]; v1 = jJet[u1, u2];
   w0 = rho^-4 rhob^2 v0; w1 = rho^-4 rhob^2 v1;
   target = -16 Pi RR hs (-rho^8 rhob/(Sqrt[2] RR^4)
      lJet[-1, w0, w1]);
   ordering = {"inner J", "outer L_-1"},
  3,
   u0 = rho^-2 rhob f0; u1 = rho^-2 rhob f1;
   u2 = rho^-2 rhob f2;
   v0 = jJet[u0, u1]; v1 = jJet[u1, u2];
   w0 = rho^-4 v0; w1 = rho^-4 v1;
   target = -16 Pi RR hs (-rho^8 rhob/(2 RR^4) jJet[w0, w1]);
   ordering = {"inner J", "outer J"},
  4,
   u0 = rho^-2 rhob^-2 f0; u1 = rho^-2 rhob^-2 f1;
   u2 = rho^-2 rhob^-2 f2;
   v0 = lJet[-1, u0, u1]; v1 = lJet[-1, u1, u2];
   w0 = rho^-4 rhob^2 RR^2 v0; w1 = rho^-4 rhob^2 RR^2 v1;
   target = -16 Pi RR hs (-rho^8 rhob/(Sqrt[2] RR^4)
      jJet[w0, w1]);
   ordering = {"inner L_-1", "outer J"}];
 residual = generated - target;
 expandedResidual = Expand[residual];
 timeRemnants = Cases[generated,
   Derivative[kt_, __][ff][__] /; kt > 0, Infinity];
 undefinedJetRemnants = Cases[generated,
   (fJet3 | fJet4 | dotJet | ddotJet)[__], Infinity];
 ass = Element[{mass, spin, len, RR, yy, mmode}, Reals] &&
   mass > 0 && len > 0 && RR > 0 && -1 < yy < 1 && hd != 0;
 <|"Stage" -> "A5 source block",
   "Block" -> block,
   "GenerationStart" -> "corresponding A4 exact block plus compactified-Delta Gaussian extension",
   "GenerationUsedTheoryFrakLJOrH" -> False,
   "Ordering" -> ordering,
   "RawResidualLeafCount" -> LeafCount[residual],
   "CompleteSourceResidual" -> FullSimplify[expandedResidual, ass],
   "ExpandedResidualExactlyZero" -> (expandedResidual === 0),
   "DiracDeltaFree" -> FreeQ[generated, DiracDelta],
   "UnevaluatedTimeDerivativeRemnants" -> timeRemnants,
   "UndefinedJetRemnants" -> undefinedJetRemnants,
   "OnlyFormalSpatialDerivativesRemain" ->
    (timeRemnants === {} && undefinedJetRemnants === {}),
   "Assumptions" -> {"M>0,L>0,R>0,-1<y<1,hatDelta!=0",
    "common f^(0,1,2) jets supplied by passed A5a",
    "fixed widths and smooth interior Gaussian"},
   "Excluded" -> {"finite-domain normalization", "axis completion",
    "turning branch", "width convergence"}|>
]
