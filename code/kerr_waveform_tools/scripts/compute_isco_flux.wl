Needs["KerrGeodesics`"];
Needs["Teukolsky`"];

inputPath = $ScriptCommandLine[[-2]];
outputPath = $ScriptCommandLine[[-1]];
input = Import[inputPath, "RawJSON"];
If[input["lmax"] =!= 6, Print["This entry point fixes lmax=6."]; Exit[2]];

chi = N[input["chi"], 17];
rISCO = N[input["r_isco"], 17];
orbit = KerrGeodesics`KerrGeoOrbit`KerrGeoOrbit[chi, rISCO, 0, 1];
fluxByL = Table[
  2 Total@Table[
    Total@Values[Teukolsky`TeukolskyMode`TeukolskyPointParticleMode[-2, ell, mm, 0, 0, orbit]["Fluxes"]["Energy"]],
    {mm, 1, ell}],
  {ell, 2, 6}];
flux = Re@Chop@Total[fluxByL];
Export[outputPath, <|
  "chi" -> chi,
  "lmax" -> 6,
  "flux_total_over_mu2" -> N[flux, 17],
  "relative_lmax_increment" -> N[Abs[Last[fluxByL]/flux], 17],
  "flux_version" -> input["flux_version"],
  "wolfram_version" -> $Version,
  "provenance" -> StringJoin["KerrGeodesics/Teukolsky compute_isco_flux.wl under ", $Version]
|>, "RawJSON"];
