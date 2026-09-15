bhptDependencyRoot="\\\\wsl.localhost\\Ubuntu\\home\\ljq\\code\\kerr-hyperboloidal\\outputs\\paper_metric_reference\\wolfram_dependencies";
Scan[PacletDirectoryLoad[FileNameJoin[{bhptDependencyRoot,#}]]&, {"SpinWeightedSpheroidalHarmonics","KerrGeodesics","Teukolsky"}];
Needs["SpinWeightedSpheroidalHarmonics`"];Needs["KerrGeodesics`"];Needs["Teukolsky`"];
FileNameJoin[{bhptDependencyRoot,#,"Kernel"}]& /@ {"SpinWeightedSpheroidalHarmonics","KerrGeodesics","Teukolsky"}
