"""Read selected radii from the original alpha=.3 Kerr or Schwarzschild flux curves in Fig.2."""
import argparse
from bisect import bisect_right
import hashlib
import json
from pathlib import Path
import pdfplumber


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pdf',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--radii',type=float,nargs='+',default=[10.,20.,30.])
    parser.add_argument('--background',choices=('kerr','schwarzschild'),default='kerr')
    args=parser.parse_args()
    color={'kerr':(.88072,.61104,.14205),'schwarzschild':(.36841,.50677,.70979)}[args.background]
    with pdfplumber.open(args.pdf) as document:
        page=document.pages[0]
        xgrid=sorted(l['x0'] for l in page.lines if l['x0']>230 and l['height']>120
                     and l['width']==0 and l['stroking_color']==.75)
        ygrid=sorted(l['top'] for l in page.lines if l['x0']>230 and l['width']>180
                     and l['height']==0 and l['stroking_color']==.75)
        curves=[c for c in page.curves if c['x0']>230 and len(c['pts'])==200
                and c['stroking_color']==color]
    if len(xgrid)!=5 or len(ygrid)!=4 or len(curves)!=2:
        raise ValueError('Unexpected original Fig.2 vector structure')
    dx=(xgrid[-1]-xgrid[0])/4;dy=(ygrid[-1]-ygrid[0])/3
    if max(abs(x-(xgrid[0]+i*dx)) for i,x in enumerate(xgrid))>1e-4:
        raise ValueError('Unexpected radial grid spacing')
    if max(abs(y-(ygrid[0]+i*dy)) for i,y in enumerate(ygrid))>1e-4:
        raise ValueError('Unexpected logarithmic flux grid spacing')
    rows=[]
    for radius in args.radii:
        x=xgrid[0]+(radius-10)*dx/10;flux={}
        for curve in curves:
            points=sorted(curve['pts'])
            if not points[0][0]<x<points[-1][0]:raise ValueError('Requested radius outside interpolable curve')
            i=bisect_right([p[0] for p in points],x)
            x0,y0=points[i-1];x1,y1=points[i]
            y=y0+(y1-y0)*(x-x0)/(x1-x0)
            branch='horizon_magnitude' if curve['dash'] else 'infinity'
            if branch in flux:raise ValueError('Ambiguous curve styles')
            value=10**(-1-(y-ygrid[0])/dy)
            flux[branch]=dict(pdf_top=y,plotted_flux=value,
                per_q2_cloud_mass=value*.3**6,bracketing_pdf_points=[points[i-1],points[i]])
        rows.append(dict(rp=radius,pdf_x=x,flux=flux))
    report=dict(status='digitized_reference_not_author_numerical_data',alpha=.3,background=args.background,
        source='https://arxiv.org/src/2501.09806v1',sha256=hashlib.sha256(args.pdf.read_bytes()).hexdigest(),
        xgrid=xgrid,xgrid_radii=[10,20,30,40,50],ygrid=ygrid,
        ygrid_fluxes=[.1,.01,.001,.0001],rows=rows,
        limitations=['Linear interpolation in PDF coordinates; no formal digitization error bound',
            'Horizon curve is a magnitude; no sign inferred from line position',
            'Conversion follows epsilon^2=alpha^6 Mc/M and the Fig.2 flux definition',
            'No solver outputs enter this reference extraction'])
    args.output.write_bytes((json.dumps(report,indent=2)+'\n').encode())
    print(json.dumps(rows))


if __name__=='__main__':main()
