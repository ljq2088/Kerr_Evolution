(* A6 targeted summary only; no heavy algebra and no full-endpoint claim. *)
Module[{base, first, chars, axis, flags},
 base = ExpandFileName@FileNameJoin[{
   DirectoryName[$InputFileName], "..", "..", "..", "..",
   "data", "kerr_point_particle_verification", "results", "a6"}];
 first = Import[FileNameJoin[{base, "a6_first_order.json"}], "RawJSON"];
 chars = Import[FileNameJoin[{base, "a6_characteristics.json"}], "RawJSON"];
 axis = Import[FileNameJoin[{base, "a6_axis_projection.json"}], "RawJSON"];
 flags = <|
   "FirstOrderIdentityChecked" -> TrueQ[first["FirstOrderIdentityChecked"]],
   "CharacteristicsChecked" -> TrueQ[chars["CharacteristicsChecked"]],
   "AxisProjectionChecked" -> TrueQ[axis["AxisProjectionChecked"]]|>;
 <|"Stage" -> "A6 targeted summary",
   "Checks" -> flags,
   "AllTargetedChecksPassed" -> And @@ Values[flags],
   "Assumptions" -> <|
    "FirstOrder" -> first["assumptions"],
    "Characteristics" -> chars["assumptions"],
    "AxisProjection" -> axis["assumptions"]|>,
   "Scope" -> {"A-F extraction and P/Q identity",
    "continuum characteristic roots and boundary directions",
    "axis indicial/Laurent, toy factor-aware basis, analytic projection"},
   "EngineeringOpen" -> {"post-ISCO transition prescription",
    "production nodes and stencil", "CFL and stability",
    "endpoint finite-limit rows", "production quadrature",
    "convergence and production readiness"},
   "FullEndpointClaimMade" -> False|>
]
