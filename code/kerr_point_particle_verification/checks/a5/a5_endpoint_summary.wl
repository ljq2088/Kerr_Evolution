(* A5 summary: imports saved staged results; no heavy algebra. *)
Module[{base, common, blocks, flags, blockResiduals, totalResidual, clean},
 base = ExpandFileName@FileNameJoin[{
   DirectoryName[$InputFileName], "..", "..", "..", "..",
   "data", "kerr_point_particle_verification", "results", "a5"}];
 common = Import[FileNameJoin[{base, "a5a_common_time_jets.json"}], "RawJSON"];
 blocks = Import[FileNameJoin[{base, #}], "RawJSON"] & /@ {
   "a5b_block1_ll.json", "a5c_block2_lj.json",
   "a5d_block3_jj.json", "a5e_block4_jl.json"};
 blockResiduals = (#1["CompleteSourceResidual"] &) /@ blocks;
 totalResidual = Total[blockResiduals];
 clean = And @@ Flatten[{
    (#1["DiracDeltaFree"] &) /@ blocks,
    (#1["OnlyFormalSpatialDerivativesRemain"] &) /@ blocks,
    ((#1["UnevaluatedTimeDerivativeRemnants"] === {}) &) /@ blocks,
    ((#1["UndefinedJetRemnants"] === {}) &) /@ blocks}];
 flags = <|
   "CommonJetsPassed" -> TrueQ[common["stage_passed"]] &&
     TrueQ[common["residuals"]["AllRequiredResidualsZero"]],
   "ExtensionWorldlineRecoveryPassed" -> TrueQ[
     common["residuals"]["Worldline"]["ExtensionRecovery"] === {0, 0, 0, 0}],
   "FourBlockEndpointsPassed" -> And @@ (# === 0 & /@ blockResiduals),
   "AllGeneratedBlocksClean" -> TrueQ[clean],
   "TotalSourceResidualZero" -> TrueQ[totalResidual === 0]|>;
 <|"Stage" -> "A5 endpoint summary",
   "RequiredFlags" -> flags,
   "BlockResiduals" -> AssociationThread[
    {"S1_LL", "S2_LJ", "S3_JJ", "S4_JL"}, blockResiduals],
   "TotalSourceResidual" -> totalResidual,
   "InteriorTimeAlgebraVerified" -> And @@ Values[flags],
   "A5FullEndpointVerifiedWithinInteriorScope" -> And @@ Values[flags],
   "GenerationUsedTheoryFrakHelpers" -> False,
   "GaussianFullLineNormalizationVerified" -> TrueQ[
    common["residuals"]["FullLineNormalization"] === <|"GR" -> 0, "Gy" -> 0|>],
   "FiniteDomainNormalizationVerified" -> False,
   "AxisCompletionVerified" -> False,
   "TurningPointBranchVerified" -> False,
   "WidthConvergenceVerified" -> False,
   "A1BoundaryInherited" -> True,
   "A4EndpointInherited" -> True|>
]
