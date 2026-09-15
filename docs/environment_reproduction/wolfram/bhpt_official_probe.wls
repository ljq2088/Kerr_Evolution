ClearAll["Global`*"];
Get[FileNameJoin[{DirectoryName[$InputFileName],"bhpt_paths.wl"}]];
outdir=FileNameJoin[{DirectoryName[$InputFileName],"bhpt_official"}];If[!DirectoryQ[outdir],CreateDirectory[outdir]];
json[x_Association]:=Association@KeyValueMap[(#1->json[#2])&,x];
json[x_List]:=json /@ x;
json[x_Complex]:=<|"re"->N[Re[x]],"im"->N[Im[x]]|>;
json[x_?NumberQ]:=N[x];
json[x_String]:=x;json[x_]:=ToString[x,InputForm];
wp=40;
parameterSets={<|"id"->"sam_a06_r8","a"->N[3/5,wp],"p"->N[8,wp],"l"->2,"spins"->{-2,-1,0,1,2},"ms"->{1},"radii"->N[{3,8,20},wp]|>,<|"id"->"actual_r20","a"->N[8771530275949366/10^16,wp],"p"->N[20,wp],"l"->1,"spins"->{-1,0,1},"ms"->{1,-1},"radii"->N[{3,20,40},wp]|>};
summary=<|"version"->$Version,"working_precision"->wp,"packages"->AssociationMap[PacletObject[#]["Version"]&,{"Teukolsky","SpinWeightedSpheroidalHarmonics","KerrGeodesics"}],"cases"->{}|>;
Export[FileNameJoin[{outdir,"progress.json"}],<|"stage"->"started"|>,"RawJSON"];
Do[
 a=param["a"];p=param["p"];ell=param["l"];radii=param["radii"];
 orbit=KerrGeoOrbit[a,p,0,1];freq=KerrGeoFrequencies[a,p,0,1];omegaPhi=Last[Values[freq]];
 orbitRecord=<|"id"->param["id"],"a"->a,"p"->p,"energy"->orbit["Energy"],"angular_momentum"->orbit["AngularMomentum"],"frequencies"->freq,"fourvelocity_at_zero"->Map[#[0]&,orbit["FourVelocity"]],"omega_phi"->omegaPhi,"raw_input_form"->ToString[{orbit["Energy"],orbit["AngularMomentum"],freq,orbit["FourVelocity"]},InputForm]|>;
 Export[FileNameJoin[{outdir,param["id"]<>"_orbit.json"}],json[orbitRecord],"RawJSON"];
 Do[
 omega=mm omegaPhi;caseid=param["id"]<>"_s"<>ToString[ss]<>"_m"<>ToString[mm];
 Export[FileNameJoin[{outdir,"progress.json"}],<|"stage"->"case_started","case"->caseid|>,"RawJSON"];
 start=AbsoluteTime[];
 result=TimeConstrained[
  sw=SpinWeightedSpheroidalHarmonicS[ss,ell,mm,a omega,Method->{"SphericalExpansion","NumTerms"->30}];
  eig=SpinWeightedSpheroidalEigenvalue[ss,ell,mm,a omega];
  angles=N[{Pi/4,Pi/2,3 Pi/4},wp];
  angular=Table[<|"theta"->th,"S"->sw[th,0],"dS"->Derivative[1,0][sw][th,0]|>,{th,angles}];
  rad=TeukolskyRadial[ss,ell,mm,a,omega,Method->"MST",WorkingPrecision->wp,AccuracyGoal->25,PrecisionGoal->25];
  radial=Association@Table[rr=rad[bc];bc-><|"method"->rr["Method"],"eigenvalue"->rr["Eigenvalue"],"nu"->rr["RenormalizedAngularMomentum"],"amplitudes"->rr["Amplitudes"],"samples"->Table[<|"r"->rv,"R"->rr[rv],"dR"->rr'[rv]|>,{rv,radii}]|>,{bc,{"In","Up"}}];
  raw=<|"case"->caseid,"a"->a,"p"->p,"s"->ss,"l"->ell,"m"->mm,"omega"->omega,"eigenvalue"->eig,"angular"->angular,"radial"->radial|>;
  <|"ok"->True,"data"->json[raw],"raw_input_form"->ToString[raw,InputForm]|>,120,<|"ok"->False,"reason"->"TimeConstrained120s"|>];
 result=Join[result,<|"elapsed_seconds"->N[AbsoluteTime[]-start]|>];
 Export[FileNameJoin[{outdir,caseid<>".json"}],result,"RawJSON"];
 summary["cases"]=Append[summary["cases"],<|"case"->caseid,"ok"->result["ok"],"elapsed_seconds"->result["elapsed_seconds"]|>];
 Export[FileNameJoin[{outdir,"summary.json"}],summary,"RawJSON"];
 ,{ss,param["spins"]},{mm,param["ms"]}];
,{param,parameterSets}];
Export[FileNameJoin[{outdir,"progress.json"}],<|"stage"->"complete"|>,"RawJSON"];
Exit[];
