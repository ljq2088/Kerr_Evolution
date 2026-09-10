"""Read v1 Fig. 6 vector markers; requires optional pdfplumber.

Axis calibration is from visual inspection of the original complete page.
These are digitized graphic positions, not author-supplied numerical data.
"""
import argparse
import hashlib
import json
from pathlib import Path


def main():
    import pdfplumber
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pdf',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    with pdfplumber.open(args.pdf) as document:
        page=document.pages[0]
        # Original PDF coordinates in points; top measured downwards.
        rows=[]
        for curve in page.curves:
            color=curve.get('non_stroking_color')
            if not isinstance(color,(tuple,list)) or len(color)!=3:
                continue
            x=(curve['x0']+curve['x1'])/2
            y=(curve['top']+curve['bottom'])/2
            ell=round(2+(x-61.937075)/15.398085)
            if not 2<=ell<=12:
                continue
            flux=10**(-8-4*(y-38.66878)/(79.39669-38.66878))
            rows.append(dict(ell=ell,pdf_marker_center=[x,y],plotted_flux=flux))
        for ell in range(2,13):
            ordered=sorted((r for r in rows if r['ell']==ell),key=lambda r:-r['plotted_flux'])
            for rank,row in enumerate(ordered):
                row['m']=ell-2*rank
    report=dict(source='https://arxiv.org/src/2501.09806v1',
                figure='Figures/Flux_l_mode_convergence_inf.pdf',
                sha256=hashlib.sha256(args.pdf.read_bytes()).hexdigest(),
                status='digitized_vector_graphic_not_author_numeric_data',
                parameters=dict(alpha=.3,rp_over_M=20),
                mode_assignment='Largest marker m=ell, then m decreases by two; paper Fig.6 caption.',
                normalization='Axis labeled F^{s,infinity}_{ell m}; comparison to epsilon definition unresolved.',
                markers=sorted(rows,key=lambda r:(r['ell'],r['m'])))
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps([r for r in report['markers'] if r['ell']==3]))


if __name__=='__main__':
    main()
