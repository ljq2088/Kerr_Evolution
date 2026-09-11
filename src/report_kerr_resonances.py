"""Compare Kerr bound-state poles to forced-channel frequencies, without fitting flux."""
import json
from pathlib import Path
from environment_source import ThresholdCloud
from environment_bound_spectrum import quasibound_mode,continue_fundamental


def main():
    cloud=ThresholdCloud(alpha=.3);cases=[]
    for ell,m in ((0,0),(2,2)):
        history=continue_fundamental(cloud.a,.3,ell,m)
        coarse=history[-1];seed=complex(*coarse['omega'])
        refined=quasibound_mode(cloud.a,.3,ell,m,seed=seed,outer_efolds=60.,offset=1e-6,rtol=2e-12)
        omega=complex(*refined['omega'])
        op=(omega.real-cloud.omega)/(m-1)
        radius=(1/op-cloud.a)**(2/3) if 0<op<1/cloud.a else None
        drive=cloud.omega+(m-1)/(20**1.5+cloud.a)
        case=dict(ell=ell,m=m,continuation=history,refined=refined,
            frequency_change=abs(omega-seed),resonance_radius_from_real_frequency=radius,
            drive_at_rp20=drive,detuning_at_rp20=drive-omega.real,
            detuning_over_absolute_imaginary_part=(drive-omega.real)/abs(omega.imag))
        cases.append(case);print({k:v for k,v in case.items() if k not in ('continuation','refined')},omega,flush=True)
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/kerr_resonance_poles.json'
    out.write_text(json.dumps(dict(status='bound_pole_diagnostic_not_forced_flux_reproduction',
        a=cloud.a,alpha=.3,cloud_frequency=cloud.omega,cases=cases,
        limitation='A real-frequency crossing is not a prediction of the flux maximum or a floating orbit. Source overlaps, finite widths, growth, and depletion still matter.'),indent=2)+'\n')


if __name__=='__main__':main()
