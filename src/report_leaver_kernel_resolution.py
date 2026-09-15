"""Audit cancellation-sensitive Up/W series convergence at large radii."""
import json
from pathlib import Path
import mpmath as mp
from paper_leaver_radial import LeaverIngoing
from report_leaver_green_kernel import decaying_logder
from environment_radial import RadialGreen


def main():
    root=Path(__file__).resolve().parents[1];folder=root/'docs/environment_reproduction'
    old=json.loads((folder/'fresh_20260915_L18_mg-1_sl0.json').read_text());p=old['parameters']
    a,mu,w=p['metric']['a'],p['alpha'],p['omega'];far=old['samples'][-1]['r']
    full=LeaverIngoing(a,mu,w,0,0,terms=120000,dps=90);coeff=full.coeff
    green=RadialGreen(a,mu,w,0,0,rmax=1000,offset=1e-4,rtol=1e-11)
    rows=[]
    with mp.workdps(90):
        aa,mm,ww=[mp.mpf(str(v)) for v in (a,mu,w)];match=mp.mpf(250)
        logder,last=decaying_logder(match,aa,mm,ww,full.angular+aa*aa*ww*ww,40)
        for n in (60000,90000,120000):
            full.coeff=coeff[:n+1];inside,derivative=full.state_mp(match)
            amps=mp.lu_solve(mp.matrix([[inside,mp.conj(inside)],[derivative,mp.conj(derivative)]]),mp.matrix([1,logder]))
            W=(match-full.rp)*(match-full.rm)*(inside*logder-derivative)
            ri,_=full.state_mp(far);kernel=(amps[0]*ri+amps[1]*mp.conj(ri))/W
            reference=green.upsol.sol(far)[0]/green.w0
            row=dict(terms=n,r=far,kernel=[mp.nstr(mp.re(kernel),40),mp.nstr(mp.im(kernel),40)],
                relative_to_production=float(abs(complex(kernel)/reference-1)))
            rows.append(row);print(json.dumps(row),flush=True)
    with mp.workdps(60):
        fine=mp.mpc(*rows[-1]['kernel'])
        for row in rows:
            z=mp.mpc(*row['kernel']);row['relative_to_120000']=float(abs(z/fine-1))
    if rows[-1]['relative_to_production']>1e-8:raise RuntimeError('Full series unresolved at farthest source node')
    result=dict(status='Up_over_W_far_node_series_resolution',a=a,mu=mu,omega=w,match=250,asymptotic_order=40,dps=90,rows=rows,
        explanation='Up is obtained by cancellation of growing In and conjugate In. Relative accuracy of In alone does not guarantee accuracy of Up/W.',
        limitations=['Checks largest source radius for one bound scalar00 mode; no all-mode coverage',
          'The insufficient 60000-term diagnostic never generated historical production fluxes'])
    (folder/'leaver_kernel_resolution_audit.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
