(* W02.2/A4: fixed-m, y=-Cos[theta], and exact distributional source. *)
Module[{kinematics, sources, group},
 group = If[ValueQ[a4Group], a4Group, "All"];
 kinematics = If[MemberQ[{"All", "Kinematics"}, group], Module[
   {ass, dy, mode, project, azTh, angTh, angY, lTh, lY, jTh, jY,
    oTh, oY, fieldRes, angRes, lRes, jRes, dPhiRes, der1Res,
    der2Res, measureRes, deltaJacRes, weightRes, jAng, uOld, vOld,
    uNew, vNew, contractRes, tensorRes, fourierMatrix,
    fourierCoeffRes, modeCommRes, nonzeroIntegerFourierRes,
    zeroIntegerFourierRes, d22, d44, y22, y44, zeroOrder,
    poleOrders, polePowerResiduals},
   ass = Element[{mass, spin, len, rad, yy}, Reals] && mass > 0 &&
     len > 0 && rad > 0 && -1 < yy < 1;
   dy = 1 - yy^2;
   mode = q[TT, rad, -Cos[theta]] Exp[I mmode phi];
   project[z_] := FullSimplify[
     (z/Exp[I mmode phi]) /. theta -> ArcCos[-yy],
     ass && Element[{mmode, sw}, Integers]];

   azTh[z_] := -I D[z, phi] + sw Cos[theta] z;
   angTh[z_] := Csc[theta] D[Sin[theta] D[z, theta], theta] + sw z -
     Csc[theta]^2 azTh[azTh[z]];
   angY[z_] := dy D[z, {yy, 2}] - 2 yy D[z, yy] -
     (mmode - sw yy)^2/dy z + sw z;
   lTh[z_, ss_] := D[z, theta] - I Csc[theta] D[z, phi] -
     I spin Sin[theta] D[z, TT] + ss Cot[theta] z;
   lY[z_, ss_] := Sqrt[dy] D[z, yy] -
     I spin Sqrt[dy] D[z, TT] + (mmode - ss yy)/Sqrt[dy] z;
   jTh[z_] := -(2 + 4 mass rad/len^2) D[z, TT] -
     (rad^2/len^2) D[z, rad];
   jY[z_] := -(2 + 4 mass rad/len^2) D[z, TT] -
     (rad^2/len^2) D[z, rad];

   oTh[z_] :=
     (8 mass (2 mass - spin^2 rad/len^2) (1 + 2 mass rad/len^2) -
       spin^2 Sin[theta]^2) D[z, {TT, 2}]
     - 2 (len^2 - (8 mass^2 - spin^2) rad^2/len^2 +
       4 spin^2 mass rad^3/len^4) D[z, TT, rad]
     - (len^2 - 2 mass rad + spin^2 rad^2/len^2) rad^2/len^2
       D[z, {rad, 2}] - angTh[z]
     + 2 spin (1 + 4 mass rad/len^2) D[z, TT, phi]
     + 2 spin rad^2/len^2 D[z, rad, phi]
     + 2 (2 mass (-sw + 2 (sw + 2) mass rad/len^2 -
       3 spin^2 rad^2/len^4) - spin^2 rad/len^2 +
       I sw spin Cos[theta]) D[z, TT]
     + 2 rad (-(1 + sw) + (sw + 3) mass rad/len^2 -
       2 spin^2 rad^2/len^4) D[z, rad]
     + 2 spin rad/len^2 D[z, phi]
     + 2 ((1 + sw) mass rad/len^2 - spin^2 rad^2/len^4) z;

   oY[z_] :=
     (8 mass (2 mass - spin^2 rad/len^2) (1 + 2 mass rad/len^2) -
       spin^2 dy) D[z, {TT, 2}]
     - 2 (len^2 - (8 mass^2 - spin^2) rad^2/len^2 +
       4 spin^2 mass rad^3/len^4) D[z, TT, rad]
     - (len^2 - 2 mass rad + spin^2 rad^2/len^2) rad^2/len^2
       D[z, {rad, 2}] - angY[z]
     + 2 I spin mmode (1 + 4 mass rad/len^2) D[z, TT]
     + 2 I spin mmode rad^2/len^2 D[z, rad]
     + 2 (2 mass (-sw + 2 (sw + 2) mass rad/len^2 -
       3 spin^2 rad^2/len^4) - spin^2 rad/len^2 -
       I sw spin yy) D[z, TT]
     + 2 rad (-(1 + sw) + (sw + 3) mass rad/len^2 -
       2 spin^2 rad^2/len^4) D[z, rad]
     + 2 I spin mmode rad/len^2 z
     + 2 ((1 + sw) mass rad/len^2 - spin^2 rad^2/len^4) z;

   dPhiRes = FullSimplify[
     project[D[mode, phi]] - I mmode q[TT, rad, yy],
     ass && Element[{mmode, sw}, Integers]];
   angRes = FullSimplify[Together[
     project[angTh[mode]] - angY[q[TT, rad, yy]]],
     ass && Element[{mmode, sw}, Integers]];
   lRes = FullSimplify[Together[
     project[lTh[mode, ss]] - lY[q[TT, rad, yy], ss]],
     ass && Element[{mmode, sw, ss}, Integers]];
   jRes = FullSimplify[Together[
     project[jTh[mode]] - jY[q[TT, rad, yy]]],
     ass && Element[{mmode, sw}, Integers]];
   fieldRes = FullSimplify[Together[
     project[oTh[mode]] - oY[q[TT, rad, yy]]],
     ass && Element[{mmode, sw}, Integers]];

   der1Res = FullSimplify[
     (D[f[-Cos[theta]], theta] /. theta -> ArcCos[-yy]) -
       Sqrt[dy] f'[yy], ass];
   der2Res = FullSimplify[
     (Csc[theta] D[Sin[theta] D[f[-Cos[theta]], theta], theta] /.
       theta -> ArcCos[-yy]) - (dy f''[yy] - 2 yy f'[yy]), ass];
   measureRes = FullSimplify[
     (sigp rp0^2 Sin[thp]/len^2) (1/Sin[thp]) -
       sigp rp0^2/len^2,
     Element[{sigp, rp0, len, thp}, Reals] && len > 0 && rp0 > 0 &&
       sigp > 0 && 0 < thp < Pi];
   deltaJacRes = FullSimplify[Sin[thp] (1/Sin[thp]) - 1,
     0 < thp < Pi];
   weightRes = FullSimplify[
     mu len^2/(sigp rp0^2 Sin[thp] uT) Sin[thp]
       Exp[-I mmode phip]/(2 Pi) -
     mu len^2 Exp[-I mmode phip]/(2 Pi sigp rp0^2 uT),
     Element[{mu, len, sigp, rp0, thp, uT, phip}, Reals] &&
       len > 0 && sigp > 0 && rp0 > 0 && uT != 0 && 0 < thp < Pi];

   jAng = DiagonalMatrix[{1, 1, Sin[thp], 1}];
   uOld = {uT, uR, uTh, uP};
   vOld = {vT, vR, vTh, vP};
   uNew = jAng.uOld;
   vNew = Transpose[Inverse[jAng]].vOld;
   contractRes = FullSimplify[vNew.uNew - vOld.uOld, 0 < thp < Pi];
   tensorRes = FullSimplify[
     jAng.Outer[Times, uOld, uOld].Transpose[jAng] -
       Outer[Times, uNew, uNew], 0 < thp < Pi];

   fourierMatrix = Table[FullSimplify[
     Integrate[Exp[I (n - m) phi], {phi, 0, 2 Pi}]/(2 Pi) -
       KroneckerDelta[n, m]], {m, -6, 6}, {n, -6, 6}];
   fourierCoeffRes = Table[FullSimplify[
     (1/(2 Pi)) Sum[Exp[-I n phip]
       Integrate[Exp[I (n - m) phi], {phi, 0, 2 Pi}]/(2 Pi),
       {n, -8, 8}] - Exp[-I m phip]/(2 Pi)], {m, {2, 4}}];
   nonzeroIntegerFourierRes = FullSimplify[
     (Exp[2 Pi I kk] - 1)/(2 Pi I kk),
     Element[kk, Integers] && kk != 0];
   zeroIntegerFourierRes = Integrate[1, {phi, 0, 2 Pi}]/(2 Pi) - 1;
   modeCommRes = {
     FullSimplify[project[D[b[TT, rad, theta] mode, TT]] -
       D[b[TT, rad, ArcCos[-yy]] q[TT, rad, yy], TT],
       ass && Element[{mmode, sw}, Integers]],
     FullSimplify[project[D[b[TT, rad, theta] mode, rad]] -
       D[b[TT, rad, ArcCos[-yy]] q[TT, rad, yy], rad],
       ass && Element[{mmode, sw}, Integers]]};

   d22 = FullSimplify[WignerD[{2, 2, 2}, theta], 0 < theta < Pi];
   d44 = FullSimplify[WignerD[{4, 4, 2}, theta], 0 < theta < Pi];
   y22 = FullSimplify[d22 /. theta -> ArcCos[-yy], -1 < yy < 1];
   y44 = FullSimplify[d44 /. theta -> ArcCos[-yy], -1 < yy < 1];
   zeroOrder[expr_, point_] := First@Select[Range[0, 8],
     FullSimplify[Limit[D[expr, {yy, #}], yy -> point,
       Direction -> If[point == -1, "FromAbove", "FromBelow"]]] =!= 0 &];
   poleOrders = {{zeroOrder[y22, -1], zeroOrder[y22, 1]},
     {zeroOrder[y44, -1], zeroOrder[y44, 1]}};
   polePowerResiduals = poleOrders - {{0, 2}, {1, 3}};

   <|"IndependentTransformedTargets" -> <|
      "DPhiToIMResidual" -> dPhiRes,
      "AngularOperatorResidual" -> angRes,
      "SourceLResidual" -> lRes,
      "SourceJResidual" -> jRes,
      "CompleteFieldOperatorResidual" -> fieldRes|>,
     "AuxiliaryChecks" -> <|
      "FourierNonzeroIntegerResidual" -> nonzeroIntegerFourierRes,
      "FourierZeroIntegerResidual" -> zeroIntegerFourierRes,
      "FourierOrthogonalityRangeMinus6To6AllZero" ->
       And @@ (# === 0 & /@ Flatten[fourierMatrix]),
      "PointDeltaModeCoefficientResidualsM2M4" -> fourierCoeffRes,
      "ModeCommutationResiduals" -> modeCommRes,
      "FirstAngularDerivativeResidual" -> der1Res,
      "AngularDivergenceResidual" -> der2Res,
      "MeasureResidual" -> measureRes,
      "AngularDeltaJacobianResidual" -> deltaJacRes,
      "FixedMProjectionWeightResidual" -> weightRes,
      "VectorCovectorContractionResidual" -> contractRes,
      "RankTwoTensorJacobianResiduals" -> tensorRes|>,
     "ExplicitWignerDHarmonics" -> <|
      "d22" -> d22, "d44" -> d44,
      "y22" -> Factor[y22], "y44" -> Factor[y44],
      "PoleOrdersNorthSouthM2M4" -> poleOrders,
      "PolePowerResiduals" -> polePowerResiduals|>|>
   ], Missing["NotExecuted", "Kinematics"]];

 sources = If[MemberQ[{"All", "Sources"}, group], Module[
   {ass, dy, rr, del, rhoTh, rhobTh, rhoY, rhobY, project, lTh, lY,
    jTh, jY, nnTh, mnTh, mmTh, nnY, mnY, mmY, bTh, bY, blockRes,
    totalRes, op, tr, lagrange, lCoef, jCoef, lTransposeRes,
    jTransposeRes, weakDerivativeRes, forwardOrders, transposeOrders,
    expectedTransposeOrders, orderResiduals},
   dy = 1 - yy^2;
   rr = len^2/rad;
   del = rr^2 - 2 mass rr + spin^2;
   ass = Element[{mass, spin, len, rad, yy}, Reals] && mass > 0 &&
     len > 0 && rad > 0 && -1 < yy < 1 && Element[mmode, Integers];
   rhoTh = 1/(rr - I spin Cos[theta]);
   rhobTh = 1/(rr + I spin Cos[theta]);
   rhoY = 1/(rr + I spin yy);
   rhobY = 1/(rr - I spin yy);
   project[z_] := FullSimplify[
     (z/Exp[I mmode phi]) /. theta -> ArcCos[-yy], ass];
   lTh[z_, ss_] := D[z, theta] - I Csc[theta] D[z, phi] -
     I spin Sin[theta] D[z, TT] + ss Cot[theta] z;
   lY[z_, ss_] := Sqrt[dy] D[z, yy] -
     I spin Sqrt[dy] D[z, TT] + (mmode - ss yy)/Sqrt[dy] z;
   jTh[z_] := -(2 + 4 mass rad/len^2) D[z, TT] -
     (rad^2/len^2) D[z, rad];
   jY[z_] := -(2 + 4 mass rad/len^2) D[z, TT] -
     (rad^2/len^2) D[z, rad];

   nnTh = fNN[TT, rad, -Cos[theta]] Exp[I mmode phi];
   mnTh = fMN[TT, rad, -Cos[theta]] Exp[I mmode phi];
   mmTh = fMM[TT, rad, -Cos[theta]] Exp[I mmode phi];
   nnY = fNN[TT, rad, yy];
   mnY = fMN[TT, rad, yy];
   mmY = fMM[TT, rad, yy];
   bTh = {
     -rhoTh^8 rhobTh lTh[
       rhoTh^-4 lTh[rhoTh^-2 rhobTh^-1 nnTh, 0], -1],
     -(1/Sqrt[2]) rhoTh^8 rhobTh del^2 lTh[
       rhoTh^-4 rhobTh^2 jTh[
         rhoTh^-2 rhobTh^-2 del^-1 mnTh], -1],
     -(1/2) rhoTh^8 rhobTh del^2 jTh[
       rhoTh^-4 jTh[rhoTh^-2 rhobTh mmTh]],
     -(1/Sqrt[2]) rhoTh^8 rhobTh del^2 jTh[
       rhoTh^-4 rhobTh^2 del^-1 lTh[
         rhoTh^-2 rhobTh^-2 mnTh, -1]]};
   bY = {
     -rhoY^8 rhobY lY[rhoY^-4 lY[rhoY^-2 rhobY^-1 nnY, 0], -1],
     -(1/Sqrt[2]) rhoY^8 rhobY del^2 lY[
       rhoY^-4 rhobY^2 jY[rhoY^-2 rhobY^-2 del^-1 mnY], -1],
     -(1/2) rhoY^8 rhobY del^2 jY[
       rhoY^-4 jY[rhoY^-2 rhobY mmY]],
     -(1/Sqrt[2]) rhoY^8 rhobY del^2 jY[
       rhoY^-4 rhobY^2 del^-1 lY[rhoY^-2 rhobY^-2 mnY, -1]]};
   blockRes = MapThread[
     FullSimplify[Together[project[#1] - #2], ass] &, {bTh, bY}];
   totalRes = FullSimplify[Together[project[Total[bTh]] - Total[bY]], ass];

   op[z_, c_] := c[[1]] D[z, TT] + c[[2]] D[z, rad] +
     c[[3]] D[z, yy] + c[[4]] z;
   tr[z_, c_] := -D[c[[1]] z, TT] - D[c[[2]] z, rad] -
     D[c[[3]] z, yy] + c[[4]] z;
   lagrange[f_, chi_, c_] := FullSimplify[
     chi op[f, c] - f tr[chi, c] -
       (D[c[[1]] f chi, TT] + D[c[[2]] f chi, rad] +
         D[c[[3]] f chi, yy]), ass];
   lCoef = {-I spin Sqrt[dy], 0, Sqrt[dy],
     (mmode - ss yy)/Sqrt[dy]};
   jCoef = {-(2 + 4 mass rad/len^2), -rad^2/len^2, 0, 0};
   lTransposeRes = lagrange[ff[TT, rad, yy], chi[TT, rad, yy], lCoef];
   jTransposeRes = lagrange[ff[TT, rad, yy], chi[TT, rad, yy], jCoef];
   weakDerivativeRes = Table[FullSimplify[
     Nest[-D[#, rad] &, Nest[-D[#, yy] &, chi[TT, rad, yy], k], j] -
       (-1)^(j + k) D[chi[TT, rad, yy], {rad, j}, {yy, k}]],
     {j, 0, 2}, {k, 0, 2}];

   forwardOrders = {{"L0", "L-1"}, {"J", "L-1"},
     {"J-inner", "J-outer"}, {"L-1", "J"}};
   transposeOrders = Reverse /@ forwardOrders;
   expectedTransposeOrders = {{"L-1", "L0"}, {"L-1", "J"},
     {"J-outer", "J-inner"}, {"J", "L-1"}};
   orderResiduals = MapThread[SameQ, {transposeOrders, expectedTransposeOrders}];

   <|"IndependentWholeBlockTargets" -> <|
      "BlockResiduals" -> AssociationThread[
       {"nn_LL", "mbar_n_LJ", "mbar_mbar_JJ", "mbar_n_JL"}, blockRes],
      "TotalFourBlockResidual" -> totalRes,
      "AllBlockResidualsZero" -> And @@ (# === 0 & /@ blockRes)|>,
     "WeakActionAuxiliaryChecks" -> <|
      "LFormalTransposeLagrangeResidual" -> lTransposeRes,
      "JFormalTransposeLagrangeResidual" -> jTransposeRes,
      "PureDerivativeSignsJ0To2K0To2" -> weakDerivativeRes,
      "ForwardInnerToOuterOrders" -> forwardOrders,
      "TransposeActionOrders" -> transposeOrders,
      "OrderingChecks" -> orderResiduals|>|>
   ], Missing["NotExecuted", "Sources"]];

 <|"ExecutionGroup" -> group,
   "A4FixedMYDistribution" -> <|
    "KinematicsAndField" -> kinematics,
    "ExactDistributionalSources" -> sources|>,
   "InheritedOpenBoundaries" -> {
    "A1 signature/Weyl-scalar bridge",
    "A1 sourced unit normalization",
    "A1 independent Ripley-PDF transcription"},
   "GaussianAssumptionsIntroduced" -> False,
   "PowerExpandUsed" -> False|>
]
