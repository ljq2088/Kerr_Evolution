(* A4 endpoint summary: imports saved stage results; no heavy algebra. *)
Module[{base, a4a, a4b, old, flags, weakFlags, poleFlag, full},
 base = ExpandFileName@FileNameJoin[{
   DirectoryName[$InputFileName], "..", "..", "..", "..",
   "data", "kerr_point_particle_verification", "results"}];
 a4a = Import[FileNameJoin[{base, "a4", "a4a_fixed_m_theta.json"}], "RawJSON"];
 a4b = Import[FileNameJoin[{base, "a4", "a4b_y_field_source.json"}], "RawJSON"];
 old = Import[FileNameJoin[{base, "A4_fixed_m_y_distribution.json"}], "RawJSON"];
 flags = <|
   "ThetaFieldVerified" -> TrueQ[a4a["residuals"]["CompleteField"] == 0],
   "ThetaFourBlocksVerified" -> TrueQ[
     And @@ (# == 0 & /@ Values[a4a["residuals"]["SourceBlocks"]])],
   "ThetaEquationVerified" -> TrueQ[
     a4a["residuals"]["CompleteThetaEquationResidual"] == 0],
   "YFieldVerified" -> TrueQ[a4b["residuals"]["CompleteYField"] == 0],
   "YFourBlocksVerified" -> TrueQ[
     And @@ (# == 0 & /@ Values[a4b["residuals"]["CompleteSourceBlocks"]])],
   "YSingleEndToEndVerified" -> TrueQ[
     a4b["residuals"]["SingleEndToEndResidual"] == 0],
   "A4aStagePassed" -> TrueQ[a4a["stage_passed"]],
   "A4bStagePassed" -> TrueQ[a4b["stage_passed"]]|>;
 weakFlags = <|
   "LTransposeResidualZero" -> TrueQ[
     old["same_definition_auxiliary_checks"]["L_formal_transpose_lagrange"] == 0],
   "JTransposeResidualZero" -> TrueQ[
     old["same_definition_auxiliary_checks"]["J_formal_transpose_lagrange"] == 0],
   "DerivativeSignResidualsZero" -> TrueQ[
     Flatten[old["same_definition_auxiliary_checks"][
       "pure_derivative_signs_j0_to2_k0_to2"]] == ConstantArray[0, 9]],
   "BlockOrderingChecksTrue" -> TrueQ[And @@
     old["same_definition_auxiliary_checks"]["block_ordering_checks"]]|>;
 poleFlag = TrueQ[Flatten[old["explicit_wigner_d_pole_check"][
   "pole_power_residuals"]] == ConstantArray[0, 4]];
 full = And @@ Values[flags];
 <|"Stage" -> "A4 endpoint summary",
   "PrimaryEndpointFlags" -> flags,
   "A4FullEndpointVerified" -> full,
   "WeakDistributionDiagnostics" -> weakFlags,
   "PolePowerDiagnosticVerified" -> poleFlag,
   "OldArbitraryProjectionBlocksUsedAsPassEvidence" -> False,
   "HeavyAlgebraRerun" -> False,
   "Transformation" ->
    "complete periodic point-particle E0 -> fixed-m theta -> fixed-m y",
   "GaussianIntroduced" -> False,
   "SpatialDeltaJetsExpanded" -> False,
   "A1SignatureUnitPDFBoundaryInherited" -> True,
   "ProjectGateStillConditionalOnA1" -> True|>
]
