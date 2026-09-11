"""Read Fig.4 vector curves using independently checked main/inset axes."""
import argparse
import bisect
import hashlib
import json
from pathlib import Path
import pdfplumber


def sample(points,x,origin,step,frame):
    points=sorted(points)
    i=bisect.bisect_left([p[0] for p in points],x)
    if i==0 or i==len(points):raise ValueError('Requested point not bracketed')
    left,right=points[i-1:i+1]
    y=left[1]+(right[1]-left[1])*(x-left[0])/(right[0]-left[0])
    return dict(value=10**(origin[1]+(y-origin[0])*step),pdf_x=x,pdf_top=y,bracket=[left,right],
                inside_axes=frame['x0']<=x<=frame['x1'] and frame['top']<=y<=frame['bottom'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pdf',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    with pdfplumber.open(args.pdf) as pdf:
        page=pdf.pages[0]
        mainframe=next(r for r in page.rects if abs(r['x0']-44.59)<1e-5)
        insetframe=next(r for r in page.rects if abs(r['x0']-141.52602)<1e-5)
        def xticks(top):
            return sorted(l['x0'] for l in page.lines if l['x0']==l['x1'] and abs(l['top']-top)<1e-5)
        xmain=xticks(122.72488);xinset=xticks(47.69597)
        assert len(xmain)==len(xinset)==5
        for ticks in (xmain,xinset):
            assert max(abs((ticks[i]-ticks[0])-i*(ticks[-1]-ticks[0])/4) for i in range(5))<2e-5
        def yticks(x):
            return sorted(l['top'] for l in page.lines if abs(l['x0']-x)<1e-5
                          and l['top']==l['bottom'] and 3<l['x1']-l['x0']<4)
        ymain=yticks(44.59);yinset=yticks(141.52602)
        assert len(ymain)==4 and len(yinset)==3
        assert max(abs(ymain[i]-ymain[0]-i*(ymain[-1]-ymain[0])/3) for i in range(4))<3e-5
        assert abs(yinset[1]-(yinset[0]+yinset[2])/2)<1e-5
        main_step=-12/(ymain[-1]-ymain[0]);inset_step=-4/(yinset[-1]-yinset[0])
        curves={}
        for c in page.curves:
            if len(c['pts'])<900:continue
            color=c['stroking_color']
            if color==.5:kind='gravitational'
            elif tuple(color)==(.36841,.50677,.70979):kind='scalar'
            elif tuple(color)==(.92252,.38562,.20917):kind='ratio'
            else:raise ValueError('Unrecognized curve color')
            boundary='horizon_magnitude' if c.get('dash') else 'infinity'
            if (kind,boundary) in curves:raise ValueError('Duplicate curve')
            curves[(kind,boundary)]=c['pts']
        assert len(curves)==6
        rows=[]
        for radius in (10.,20.,30.,40.):
            row=dict(rp=radius)
            for kind in ('scalar','gravitational','ratio'):
                inset=kind=='ratio';xs=xinset if inset else xmain
                x=xs[0]+(radius-10)*(xs[-1]-xs[0])/40
                origin=(yinset[0],2) if inset else (ymain[0],-13)
                row[kind]={boundary:sample(curves[(kind,boundary)],x,origin,inset_step if inset else main_step,
                                          insetframe if inset else mainframe)
                           for boundary in ('infinity','horizon_magnitude')}
            rows.append(row)
    result=dict(status='digitized_reference_not_author_numerical_data',
        source='https://arxiv.org/src/2501.09806v1',pdf=args.pdf.name,
        sha256=hashlib.sha256(args.pdf.read_bytes()).hexdigest(),
        axes=dict(main_x=xmain,inset_x=xinset,main_y=ymain,inset_y=yinset),rows=rows,
        limitations=['Linear interpolation of vector coordinates in log-flux axes',
            'No rigorous digitization error bound; no numerical amplitudes used for calibration'])
    args.output.write_bytes((json.dumps(result,indent=2)+'\n').encode('utf-8'))
    print(json.dumps(rows))


if __name__=='__main__':main()
