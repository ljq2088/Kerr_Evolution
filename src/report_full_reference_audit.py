"""Independent MuPDF vector readout of both Fig.2 panels at r0=20M."""
import hashlib,json
from pathlib import Path
import numpy as np
import pymupdf as fitz


def main():
    root=Path(__file__).resolve().parents[1];folder=root/'docs/environment_reproduction'
    pdf=root/'outputs/environment_reference/Flux_Inf_Hor_TwoPanels_2.pdf'
    doc=fitz.open(pdf);drawings=doc[0].get_drawings();rows=[];input_hashes={}
    for alpha in (.2,.3):
        inpanel=lambda x:x<230 if alpha==.2 else x>230
        xs=[];ys=[]
        for d in drawings:
            if not d['color'] or max(abs(c-.75) for c in d['color'])>1e-3 or not inpanel(d['rect'].x0):continue
            for kind,p,q in d['items']:
                if kind!='l':raise ValueError('Unexpected grid primitive')
                if p.x==q.x:xs.append(p.x)
                if p.y==q.y:ys.append(p.y)
        xs=sorted(xs);ys=sorted(ys)
        if len(xs)!=5 or len(ys)!=4:raise ValueError('Unexpected grid shape')
        dx=(xs[-1]-xs[0])/4;dy=(ys[-1]-ys[0])/3
        np.testing.assert_allclose(xs,np.array([0,1,2,3,4])*dx+xs[0],atol=1e-4,rtol=0)
        for background,color in [('kerr',[.88072,.61104,.14205]),('schwarzschild',[.36841,.50677,.70979])]:
            curves=[d for d in drawings if d['color'] and max(abs(c-r) for c,r in zip(d['color'],color))<1e-5 and len(d['items'])>=20 and d['rect'].width>100 and inpanel(d['rect'].x0)]
            if len(curves)!=2:raise ValueError('Unexpected curve count')
            for d in curves:
                pts=np.array([list(d['items'][0][1])]+[list(item[2]) for item in d['items']])
                x=xs[1];y=np.interp(x,pts[:,0],pts[:,1]);plotted=10**(-1-(y-ys[0])/dy)
                branch='infinity' if d['dashes']=='[] 0' else 'horizon'
                row=dict(alpha=alpha,background=background,boundary=branch,r0=20.,plotted=plotted,per_q2_eta=plotted*alpha**6)
                suffix=('_alpha0.2' if alpha==.2 else '')+('_schwarzschild_frozen' if background=='schwarzschild' else '')
                cp=folder/f'flux_coverage_L18_nt18_i6_h5_f5{suffix}.json'
                if cp.exists():
                    coverage=json.loads(cp.read_text());block=coverage[branch];value=block['finite_resolution_total']
                    row['coverage_sha256']=hashlib.sha256(cp.read_bytes()).hexdigest()
                    if value is not None:
                        total=0.
                        for channel in block['modes']:
                            actualpath=folder/channel['file'];actual=json.loads(actualpath.read_text())
                            input_hashes[channel['file']]=hashlib.sha256(actualpath.read_bytes()).hexdigest()
                            if actual['flux']!=channel['flux']:raise ValueError('Coverage differs from actual response')
                            total+=actual['flux'][branch]['orbital_energy']
                        np.testing.assert_allclose(value,total,rtol=1e-14)
                        row.update(computed_signed=value,relative_percent=100*(abs(value)/row['per_q2_eta']-1))
                rows.append(row)
    result=dict(status='independent_pdf_readout_and_existing_finite_totals',rows=rows,
        pdf_sha256=hashlib.sha256(pdf.read_bytes()).hexdigest(),input_response_sha256=input_hashes,
        limitations=['MuPDF uses single-precision PDF coordinates; no formal figure digitization error bound',
            'All readouts are r0=20M, not an orbit sweep',
            'Alpha=.2 source cutoffs have not been separately converged',
            'Historical response files are checked, not relabeled as fresh code-versioned results'])
    (folder/'full_reference_audit_20260915.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(rows,indent=2))


if __name__=='__main__':main()
