"""Aggregate independently solved finite gravitational multipole sums."""
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    folder=Path(__file__).resolve().parents[1]/"docs/environment_reproduction"
    rows=[]
    for radius in [3.5,4,5,*range(6,51,2)]:
        path=folder/f"gravitational_flux_r{radius:g}_L12.json"
        d=json.loads(path.read_text())
        assert d["status"]=="finite_gravitational_flux_sum_not_full_paper_reproduction"
        assert d["parameters"]["rp"]==radius and d["parameters"]["ellmax"]==12
        assert d["parameters"]["a"]==.8771530275949366
        modes={(x["ell"],x["m"]):x for x in d["modes"]}
        assert len(modes)==154
        for (ell,m),v in modes.items():
            for boundary in ("infinity","horizon"):
                np.testing.assert_allclose(v[boundary],modes[ell,-m][boundary],rtol=1e-10,atol=1e-300)
        total=d["totals"]
        for boundary in total:
            np.testing.assert_allclose(total[boundary],sum(x[boundary] for x in modes.values()),rtol=1e-13)
        rows.append(dict(rp=radius,totals=total,positive_22=modes[2,2],
            ell10_to_12_relative_change={k:total[k]/d["shells"][8]["cumulative"][k]-1 for k in total},
            input_file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    result=dict(status="complete_finite_gravitational_sweep_not_full_paper_reproduction",rows=rows,
        units="per q^2",spin=.8771530275949366,
        limitations=["Finite ell<=12 sum, not an infinite multipole error bound",
            "No scalar flux curve or environmental ratio is inferred from these data"])
    (folder/"gravitational_flux_sweep_L12.json").write_text(json.dumps(result,indent=2)+"\n")
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout="constrained")
    r=[x["rp"] for x in rows]
    for ax,boundary in zip(axes,("infinity","horizon")):
        ax.semilogy(r,[abs(x["totals"][boundary]) for x in rows],"o-",ms=3,label="Sum: ell=2..12, both signs of m")
        ax.semilogy(r,[abs(x["positive_22"][boundary]) for x in rows],"--",label="Single positive (2,2)")
        ax.set(xlabel="Orbital radius r/M",ylabel="Energy flux magnitude / q²",title=boundary.capitalize())
        ax.grid(alpha=.25);ax.legend(fontsize=8)
    fig.suptitle("Vacuum Kerr gravitational flux; a/M=0.8771530")
    fig.savefig(folder/"gravitational_flux_sweep_L12.png",dpi=160)
    print(json.dumps({"radii":len(rows),"largest_ell10_to_12_change":max(abs(v) for x in rows for v in x["ell10_to_12_relative_change"].values())}))


if __name__=="__main__":main()
