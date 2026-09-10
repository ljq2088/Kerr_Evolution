(* A6 analytic continuous/discrete spin-weighted projection check. *)
Module[{swsh, sets, continuous, n, nodes, mat, rhs, weights,
  momentResidual, discrete, cres, dres},
 swsh[l_, m_, y_] := FullSimplify[
   Sqrt[(2 l + 1)/(4 Pi)] WignerD[{l, m, 2}, ArcCos[-y]], -1 < y < 1];
 sets = {{{2, 2}, {3, 2}}, {{4, 4}, {5, 4}}};
 continuous = Table[Table[FullSimplify[2 Pi Integrate[
     swsh[set[[i, 1]], set[[i, 2]], y]
      swsh[set[[j, 1]], set[[j, 2]], y], {y, -1, 1}]],
    {i, 2}, {j, 2}], {set, sets}];
 n = 14; nodes = Table[-1 + 2 j/(n + 1), {j, n}];
 mat = Table[nodes[[j]]^k, {k, 0, n - 1}, {j, n}];
 rhs = Table[If[EvenQ[k], 2/(k + 1), 0], {k, 0, n - 1}];
 weights = LinearSolve[mat, rhs];
 momentResidual = Simplify[mat.weights - rhs];
 discrete = Table[Table[FullSimplify[2 Pi Sum[
     weights[[q]] swsh[set[[i, 1]], set[[i, 2]], nodes[[q]]]
      swsh[set[[j, 1]], set[[j, 2]], nodes[[q]]], {q, n}]],
    {i, 2}, {j, 2}], {set, sets}];
 cres = Simplify[# - IdentityMatrix[2]] & /@ continuous;
 dres = Simplify[# - IdentityMatrix[2]] & /@ discrete;
 <|"Stage" -> "A6 projection normalization and leakage",
   "ContinuousGram" -> continuous, "DiscreteGram" -> discrete,
   "ContinuousResiduals" -> cres, "DiscreteResiduals" -> dres,
   "AdjacentLLeakage" -> <|"Continuous" ->
     {continuous[[1, 1, 2]], continuous[[2, 1, 2]]},
    "Discrete" -> {discrete[[1, 1, 2]], discrete[[2, 1, 2]]}|>,
   "QuadratureMomentResidual" -> momentResidual,
   "ProjectionNormalizationLeakageChecked" -> And @@ (# === 0 & /@
     Flatten[{cres, dres, momentResidual}]),
   "Assumptions" -> {"analytic Wigner-d s=-2 convention",
    "14 rational open nodes with exact interpolatory moment weights"}|>
]
