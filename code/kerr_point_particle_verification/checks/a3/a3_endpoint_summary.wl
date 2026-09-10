(* A3 endpoint summary: imports stage JSON only; no heavy algebra is rerun. *)
Module[{base, a3a, inv, bl, proj, flags, full},
 base = ExpandFileName@FileNameJoin[{
   DirectoryName[$InputFileName], "..", "..", "..", "..",
   "data", "kerr_point_particle_verification", "results", "a3"}];
 a3a = Import[FileNameJoin[{base, "a3a_geodesic_velocities.json"}], "RawJSON"];
 inv = Import[FileNameJoin[{base, "a3b_invariant_constant_t.json"}], "RawJSON"];
 bl = Import[FileNameJoin[{base, "a3b_bl_delta_pullback.json"}], "RawJSON"];
 proj = Import[FileNameJoin[{base, "a3c_tetrad_projections.json"}], "RawJSON"];
 flags = <|
   "A3aGeodesicPassed" -> TrueQ[a3a["stage_passed"]] &&
     TrueQ[a3a["residuals"]["AllRequiredResidualsZero"]],
   "InvariantConstantTRoutePassed" -> TrueQ[inv["stage_passed"]] &&
     TrueQ[inv["residuals"]["AllRequiredResidualsZero"]],
   "BLDeltaPullbackRoutePassed" -> TrueQ[bl["stage_passed"]] &&
     TrueQ[bl["residuals"]["AllRequiredResidualsZero"]],
   "TetradProjectionRoutesPassed" -> TrueQ[proj["stage_passed"]] &&
     TrueQ[proj["residuals"]["AllRequiredResidualsZero"]],
   "A2OperatorsExcluded" -> TrueQ[Not[proj["A2FourBlockOperatorsUsed"]]]|>;
 full = And @@ Values[flags];
 <|"Stage" -> "A3 endpoint summary",
   "ImportedStageFiles" -> {
    "a3a_geodesic_velocities.json",
    "a3b_invariant_constant_t.json",
    "a3b_bl_delta_pullback.json",
    "a3c_tetrad_projections.json"},
   "RequiredStageFlags" -> flags,
   "A3FullEndpointVerified" -> full,
   "A3Endpoint" ->
    {"T_nn^(H,pp)", "T_mbar_n^(H,pp)", "T_mbar_mbar^(H,pp)"},
   "A3ClosesBeforeA2FourBlockFunctional" -> True,
   "VerificationScope" -> "algebraic closure under saved local assumptions",
   "GlobalWorldlineMonotonicityVerified" -> False,
   "GlobalUniqueSliceIntersectionVerified" -> False,
   "A1SignatureTetradBridgeInherited" -> True,
   "ProjectGateStillConditionalOnA1" -> True|>
]
