"""Reconstruct Dyson Eq. (21) thresholds from independently solved cloud spectra.

This reads cached background eigenvalues, not flux or fitted wake amplitudes.
Run from any directory: .venv/bin/python src/report_paper_threshold_comparison_20260916.py
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/environment_reproduction"
STEM = "paper_threshold_comparison_20260916"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(name: str):
    return json.loads((OUT / name).read_text())


def critical_radius(a: float, omega: float, mu: float, mg: int = 1) -> float:
    return (mg / (mu - omega) - a) ** (2.0 / 3.0)


def wave_number_squared(rp, a: float, omega: float, mu: float):
    return (omega + 1.0 / (np.asarray(rp) ** 1.5 + a)) ** 2 - mu ** 2


def main() -> None:
    cloud = read("cloud_211.json")
    schw = next(x["refined"] for x in read("schwarzschild_cloud_spectrum.json")["runs"]
                if x["refined"]["alpha"] == 0.3)
    leaver_row = next(x for x in read("leaver_cloud_full_audit.json")["rows"] if x["alpha"] == 0.3)
    leaver = max(leaver_row["cases"], key=lambda x: x["terms"])
    mu = 0.3
    specifications = [
        ("Schwarzschild", 0.0, schw["omega"][0], schw["omega"][1], 41.01,
         "Independent two-sided shooting; real part of complex frequency frozen"),
        ("Kerr", cloud["refined"]["a_over_M"], cloud["refined"]["M_omega_c"], 0.0, 41.66,
         "Independent two-sided shooting of the synchronized |211> cloud"),
        ("Newtonian", 0.0, mu * (1.0 - mu ** 2 / 8.0), 0.0, 44.44,
         "Leading hydrogenic n=2 energy; not a relativistic eigenfrequency"),
    ]
    rows = []
    for name, a, omega, imag, printed, method in specifications:
        rc = critical_radius(a, omega, mu)
        row = dict(background=name, alpha=mu, cloud_n=2, cloud_ell=1, cloud_m=1,
                   scalar_m=2, metric_m=1, a_over_M=a, M_omega_cloud_real=omega,
                   M_omega_cloud_imag=imag, computed_rcrit_over_M=rc,
                   paper_printed_rcrit_over_M=printed,
                   computed_minus_printed_over_M=rc-printed,
                   paper_last_decimal_half_unit_M=0.005,
                   within_paper_rounding_bin=bool(abs(rc-printed) <= 0.005), method=method)
        assert row["within_paper_rounding_bin"]
        rows.append(row)
    kerr = rows[1]
    cases = []
    for rp in (41.6, 41.8):
        om = 1.0 / (rp ** 1.5 + kerr["a_over_M"])
        ksquared = float(wave_number_squared(rp, kerr["a_over_M"], kerr["M_omega_cloud_real"], mu))
        kval = abs(ksquared) ** 0.5
        cases.append(dict(rp_over_M=rp, M_Omega_p=om,
                          M_omega_scalar=kerr["M_omega_cloud_real"]+om,
                          M2_k_infinity_squared=ksquared,
                          asymptotic_branch="radiative" if ksquared > 0 else "evanescent",
                          M_k_infinity_real=kval if ksquared > 0 else 0.0,
                          M_k_infinity_imag=kval if ksquared < 0 else 0.0,
                          asymptotic_length_over_M=float(2*np.pi/kval) if ksquared > 0 else float(1/kval),
                          length_definition="wavelength 2*pi/k" if ksquared > 0 else "decay length 1/abs(k)"))
    assert cases[0]["M2_k_infinity_squared"] > 0 > cases[1]["M2_k_infinity_squared"]
    leaver_rc = critical_radius(float(leaver["a"]), float(leaver["omega"]), mu)
    assert abs(leaver_rc-kerr["computed_rcrit_over_M"]) < 1e-6
    sources = [OUT / "cloud_211.json", OUT / "schwarzschild_cloud_spectrum.json",
               OUT / "leaver_cloud_full_audit.json",
               ROOT / "outputs/paper_original_reference/arxiv_v1_source/main_PRL.tex"]
    report = dict(
        status="independent_background_threshold_agreement_not_full_figure5_reproduction",
        source=dict(arxiv="https://arxiv.org/abs/2501.09806v1",
                    source_url="https://arxiv.org/src/2501.09806v1",
                    tex="outputs/paper_original_reference/arxiv_v1_source/main_PRL.tex",
                    equation_tex_lines=[366,374], figure5_tex_lines=[654,662],
                    paper_equation="rcrit=[(m-mb)/(mu-Re(omega_cloud))-a]^(2/3) M"),
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in sources},
        generator_sha256=sha(Path(__file__)),
        threshold_rows=rows, figure5_orbits=cases,
        independent_cloud_crosscheck=dict(method="400-term Leaver continued fraction at 50 digits (cached)",
            a=leaver["a"],omega=leaver["omega"],continued_fraction_residual=leaver["continued_fraction_residual"],
            threshold_over_M=leaver_rc,
            threshold_difference_from_shooting_M=leaver_rc-kerr["computed_rcrit_over_M"],
            shooting_boundary_cutoff_change_M=cloud["threshold_cutoff_change"]),
        limitations=[
            "The 0.005 M band is the last printed decimal's half-unit, not a measured or rigorous uncertainty.",
            "This tests a propagation threshold, not a quasi-bound pole resonance or the flux amplitude.",
            "Schwarzschild uses Re(omega_cloud) in the paper's frozen-background threshold convention; Im(omega_cloud) is retained in metadata only.",
            "The exact Kerr spin comes from synchronization; the paper prints rounded a=0.88.",
            "No flux or field amplitude is fitted, recomputed, or claimed reproduced by this plot.",
            "The large-r lengths describe asymptotic waves; the long-range Coulomb phase is not represented by a plane wave over the entire numerical domain."
        ])
    with (OUT / f"{STEM}.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n"); writer.writeheader(); writer.writerows(rows)
    with (OUT / f"{STEM}_figure5_orbits.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(cases[0]), lineterminator="\n"); writer.writeheader(); writer.writerows(cases)
    rp_grid = np.linspace(38.0, 47.0, 1001)
    curves = {row["background"]: wave_number_squared(rp_grid, row["a_over_M"], row["M_omega_cloud_real"], mu) for row in rows}
    with (OUT / f"{STEM}_curves.csv").open("w", newline="") as f:
        writer = csv.writer(f, lineterminator="\n"); writer.writerow(["rp_over_M"]+[f"M2_k2_{name.lower()}" for name in curves])
        writer.writerows(zip(rp_grid, *curves.values()))
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(12.2, 4.9), gridspec_kw={"width_ratios":[1.28,1.0]})
    colors = {"Schwarzschild":"#3977ae", "Kerr":"#d16035", "Newtonian":"#22846f"}
    ax.axhline(0, color="#525963", lw=.9)
    ax.axhspan(0, 10, color="#eaf4ed", alpha=.5)
    ax.axhspan(-10, 0, color="#f1eef8", alpha=.5)
    for row in rows:
        name = row["background"]
        ax.plot(rp_grid, curves[name]*1e4, label=f'{name}: {row["computed_rcrit_over_M"]:.5f} M', lw=2, color=colors[name])
        ax.plot(row["computed_rcrit_over_M"],0,"o",color=colors[name],ms=4)
    for item in cases:
        rp=item["rp_over_M"]; y=item["M2_k_infinity_squared"]*1e4
        ax.plot(rp,y,"o",mfc="white",mec=colors["Kerr"],mew=1.5,ms=6,zorder=5)
        ax.annotate(f'Fig. 5: {rp:.1f} M', (rp,y), xytext=(-100,60) if rp==41.6 else (38,-65),
                    textcoords="offset points",arrowprops={"arrowstyle":"-","color":"#555555"},fontsize=9)
    ax.text(38.18,4.9,r"$k_\infty^2>0$: radiative",fontsize=9,color="#285f36")
    ax.text(38.18,-4.3,r"$k_\infty^2<0$: evanescent",fontsize=9,color="#645278")
    ax.set(xlim=(38,47),ylim=(-5.0,6.0),xlabel=r"Orbital radius $r_p/M$",ylabel=r"$10^4 M^2 k_\infty^2$",
           title=r"Scalar $m=2$ threshold ($\alpha=0.3$, cloud $m_b=1$)")
    ax.legend(loc="upper right",fontsize=8.5,frameon=True)
    bx.axvspan(-5,5,color="#f0f1f3", label="Half-unit of paper's last printed decimal")
    bx.axvline(0,color="#73777d",lw=.9,ls="--")
    for j,row in enumerate(rows):
        bx.plot((row["computed_rcrit_over_M"]-row["paper_printed_rcrit_over_M"])*1000,j,"o",ms=8,color=colors[row["background"]])
        bx.text(-5.05,j+.24, f'Computed {row["computed_rcrit_over_M"]:.8f} M  |  Paper {row["paper_printed_rcrit_over_M"]:.2f} M',fontsize=8.6)
    bx.set(yticks=range(3),yticklabels=[x["background"] for x in rows],ylim=(-.6,2.6),xlim=(-5.4,5.4),
           xlabel=r"$(r_{\rm crit}^{\rm computed}-r_{\rm crit}^{\rm printed})/(10^{-3}M)$",
           title="Agreement with the printed thresholds")
    bx.invert_yaxis(); bx.grid(axis="x",alpha=.15)
    bx.text(.5,.04,r"Gray band: $\pm0.005M$ rounding interval, not an error bar.",
            transform=bx.transAxes,ha="center",fontsize=8)
    fig.suptitle("Dyson et al.: independently reconstructed propagation thresholds",fontsize=13,y=.97)
    fig.text(.5,.025,"Eq. (21) and the mechanism of Fig. 5. Cached independent cloud spectra; no fitted field or flux amplitudes.",ha="center",fontsize=9,color="#555555")
    fig.subplots_adjust(left=.07,right=.975,top=.84,bottom=.17,wspace=.4)
    fig.savefig(OUT / f"{STEM}.png",dpi=200)
    fig.savefig(OUT / f"{STEM}.pdf")
    plt.close(fig)
    report["output_sha256"] = {str(p.relative_to(ROOT)):sha(p) for p in sorted(OUT.glob(f"{STEM}*")) if p.suffix!=".json"}
    (OUT / f"{STEM}.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"thresholds":rows,"figure5_orbits":cases,"cloud_crosscheck":report["independent_cloud_crosscheck"]},indent=2))


if __name__ == "__main__":
    main()
