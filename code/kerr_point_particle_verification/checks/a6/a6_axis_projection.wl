(* A6 axis check: indicial/Laurent behavior and factor-aware basis rows. *)
Module[
 {ang, northCoeff, southCoeff, genericResiduals, cases, laurent,
  factorResidual, northNodes, southNodes, basis, axisValues, basisValues},
 ang[z_, m_] := (1 - yy^2) D[z, {yy, 2}] - 2 yy D[z, yy] -
   (m + 2 yy)^2/(1 - yy^2) z - 2 z;
 northCoeff[p_, m_] := FullSimplify[Limit[
   xx^(1 - p) (ang[(1 + yy)^p, m] /. yy -> -1 + xx),
   xx -> 0, Direction -> "FromAbove"]];
 southCoeff[p_, m_] := FullSimplify[Limit[
   zz^(1 - p) (ang[(1 - yy)^p, m] /. yy -> 1 - zz),
   zz -> 0, Direction -> "FromAbove"]];
 genericResiduals = {
   FullSimplify[northCoeff[pp, mm] - (2 pp^2 - (mm - 2)^2/2),
    Element[{pp, mm}, Integers] && pp >= 0],
   FullSimplify[southCoeff[pp, mm] - (2 pp^2 - (mm + 2)^2/2),
    Element[{pp, mm}, Integers] && pp >= 0]};
 cases = {{2, 0, 2}, {4, 1, 3}};
 laurent = Table[With[{m = row[[1]], pn = row[[2]], ps = row[[3]]},
    Module[{ne, se, ns, ss},
     ne = ang[(1 + yy)^pn (g0 + g1 (1 + yy) + g2 (1 + yy)^2), m] /.
       yy -> -1 + xx;
     se = ang[(1 - yy)^ps (h0 + h1 (1 - yy) + h2 (1 - yy)^2), m] /.
       yy -> 1 - zz;
     ns = Normal@Series[ne, {xx, 0, 2}];
     ss = Normal@Series[se, {zz, 0, 2}];
     <|"m" -> m, "Powers" -> {pn, ps},
      "NorthIndicial" -> northCoeff[pn, m],
      "SouthIndicial" -> southCoeff[ps, m],
      "NorthNegative" -> {Coefficient[ns, xx, -2], Coefficient[ns, xx, -1]},
      "SouthNegative" -> {Coefficient[ss, zz, -2], Coefficient[ss, zz, -1]}|>
     ]], {row, cases}];
 factorResidual[nodes_, center_, p_, north_, order_] := Module[
   {vals, mat, rhs, w},
   vals = If[north, (1 + nodes)^p, (1 - nodes)^p];
   mat = Table[vals[[j]] If[k == 0, 1, (nodes[[j]] - center)^k],
     {k, 0, Length[nodes] - 1}, {j, Length[nodes]}];
   rhs = Table[D[
       If[north, (1 + x)^p, (1 - x)^p] (x - center)^k,
       {x, order}] /. x -> center, {k, 0, Length[nodes] - 1}];
   w = LinearSolve[mat, rhs];
   Simplify[mat.w - rhs]
   ];
 northNodes = {-9/10, -4/5, -7/10, -3/5, -1/2};
 southNodes = {1/2, 3/5, 7/10, 4/5, 9/10};
 basis = Table[With[{m = row[[1]], pn = row[[2]], ps = row[[3]]},
    <|"m" -> m, "NorthPower" -> pn, "SouthPower" -> ps,
     "NorthD1" -> factorResidual[northNodes, northNodes[[2]], pn, True, 1],
     "NorthD2" -> factorResidual[northNodes, northNodes[[2]], pn, True, 2],
     "SouthD1" -> factorResidual[southNodes, southNodes[[-2]], ps, False, 1],
     "SouthD2" -> factorResidual[southNodes, southNodes[[-2]], ps, False, 2]|>],
   {row, cases}];
 axisValues = Join[genericResiduals, Flatten[Table[{
     laurent[[j]]["NorthIndicial"], laurent[[j]]["SouthIndicial"],
     laurent[[j]]["NorthNegative"], laurent[[j]]["SouthNegative"]},
    {j, Length[laurent]}]]];
 basisValues = Flatten[Table[{
    basis[[j]]["NorthD1"], basis[[j]]["NorthD2"],
    basis[[j]]["SouthD1"], basis[[j]]["SouthD2"]}, {j, Length[basis]}]];
 <|"Stage" -> "A6 axis indicial and factor-aware basis",
   "GenericIndicialResiduals" -> AssociationThread[{"North", "South"},
    genericResiduals],
   "ModeLaurentChecks" -> laurent,
   "FactorAwareBasisChecks" -> basis,
   "AxisIndicialAndLaurentChecked" -> And @@ (# === 0 & /@ axisValues),
   "FactorAwareBasisChecked" -> And @@ (# === 0 & /@ basisValues),
   "Assumptions" -> {"s=-2", "m=2,4", "rational open nodes exclude endpoints"},
   "EngineeringBoundary" -> "production nodes/stencil and finite-limit rows remain W02.3"|>
]
