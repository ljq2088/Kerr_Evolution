"""Compute both Li backgrounds with 150/300-term Leaver controls before fluxes."""
import json,hashlib,time
from pathlib import Path
import mpmath as mp
from li_leaver_cloud import LiLeaverCloud
ROOT=Path(__file__).resolve().parents[1]
def pair(z):return [mp.nstr(mp.re(z),50),mp.nstr(mp.im(z),50)]
def main():
    out=ROOT/'docs/li_alignment/leaver_backgrounds.json';rows=[]
    for ell in (1,2):
        for terms in (150,300):
            start=time.perf_counter();c=LiLeaverCloud(ell=ell,m=ell,terms=terms)
            rows.append(dict(ell_c=ell,m_c=ell,n_c=0,terms=terms,dps=c.dps,
                omega=pair(c.omega),Lambda=pair(c.data['A']),cf_residual=pair(c.data['residual']),
                reduction_defect=pair(c.data['defect']),
                radial=[dict(r=r,R=pair(c.radial_mp(r)),Rprime=pair(c.radial_mp(r,1)),
                    scaled_ODE_residual=mp.nstr(c.radial_residual(r),20)) for r in (2.,3.,10.,20.,100.,200.,320.)],
                elapsed_seconds=time.perf_counter()-start))
            result=dict(status='background_spectrum_and_profile_in_progress' if len(rows)<4 else 'both_Leaver_backgrounds_computed_mass_normalization_pending',
                a='.88',alpha='.3',normalization='UNNORMALIZED a0=1; no physical mass claimed',
                frequency_treatment='Complex spectrum retained; no temporal freezing',rows=rows,
                implementation_sha256=hashlib.sha256((ROOT/'src/li_leaver_cloud.py').read_bytes()).hexdigest())
            tmp=out.with_suffix('.tmp');tmp.write_text(json.dumps(result,indent=2)+'\n');tmp.replace(out)
            print(ell,terms,rows[-1]['omega'],max(float(x['scaled_ODE_residual']) for x in rows[-1]['radial']),flush=True)
if __name__=='__main__':main()
