"""Read the eleven colored diamonds and logarithmic grid of original v1 Fig.7."""
import argparse
import hashlib
import json
import math
import pdfplumber


def main():
    from pathlib import Path
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pdf',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    with pdfplumber.open(args.pdf) as document:
        page=document.pages[0]
        diamonds=sorted([c for c in page.curves if c['fill'] and len(c['pts'])==5
                         and isinstance(c['non_stroking_color'],tuple)],key=lambda c:c['x0'])
        grid=sorted(line['top'] for line in page.lines if line['width']>180
                    and line['height']==0 and line['stroking_color']==.75)
    if len(diamonds)!=11 or len(grid)!=4:
        raise ValueError('Unexpected Fig.7 vector structure')
    step=(grid[-1]-grid[0])/3
    if max(abs(y-(grid[0]+i*step)) for i,y in enumerate(grid))>1e-4:
        raise ValueError('Nonuniform logarithmic grid')
    xs=[(d['x0']+d['x1'])/2 for d in diamonds]
    slope=(xs[-1]-xs[0])/math.log(12/2)
    if max(abs(x-(xs[0]+slope*math.log(ell/2))) for ell,x in zip(range(2,13),xs))>.01:
        raise ValueError('Markers do not match ell=2..12 on logarithmic x axis')
    rows=[dict(ell=ell,pdf_x=x,pdf_top=(d['top']+d['bottom'])/2,
        plotted_value=10**(-1-((d['top']+d['bottom'])/2-grid[0])/step))
        for ell,x,d in zip(range(2,13),xs,diamonds)]
    result=dict(status='vector_digitization_not_author_numerical_data',
        source='https://arxiv.org/src/2501.09806v1',sha256=hashlib.sha256(args.pdf.read_bytes()).hexdigest(),
        log_y_grid_top=grid,log_y_grid_values=[.1,.01,.001,.0001],markers=rows,
        limitations=['No formal digitization error bound',
            'Figure7 angular evaluation and normalization must be checked independently'])
    args.output.write_bytes((json.dumps(result,indent=2)+'\n').encode())
    print(json.dumps(rows))


if __name__=='__main__':main()
