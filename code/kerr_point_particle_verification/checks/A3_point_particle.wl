(* Compatibility dispatcher for the staged A3 verification. No heavy algebra. *)
Module[{summaryPath, summary},
 summaryPath = FileNameJoin[{
   DirectoryName[$InputFileName], "a3", "a3_endpoint_summary.wl"}];
 summary = Get[summaryPath];
 <|"CompatibilityDispatcher" -> True,
   "SupersededMonolithicEntry" -> True,
   "StagedSummaryPath" -> "checks/a3/a3_endpoint_summary.wl",
   "HeavyAlgebraRerun" -> False,
   "StagedSummary" -> summary|>
]
