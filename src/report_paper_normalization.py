"""Cross-check original v1 Fig.2 against digitized Fig.6; needs pdfplumber.

Calibration is read from the original PDF axes. No solver output enters this
comparison. Linear interpolation in PDF coordinates introduces digitization
error; the reported difference is not a numerical error bound.
"""
import argparse
import hashlib
import json
from pathlib import Path
from bisect import bisect_right


def main():
    import pdfplumber
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('figure2', type=Path)
    p.add_argument('markers', type=Path)
    p.add_argument('output', type=Path)
    args = p.parse_args()
    with pdfplumber.open(args.figure2) as pdf:
        curves = [c for c in pdf.pages[0].curves
                  if tuple(c.get('stroking_color') or ()) == (.88072, .61104, .14205)
                  and len(c['pts']) == 200 and abs(c['x0']-232.3337) < .001
                  and abs(c['x1']-417.10818) < .001]
    if len(curves) != 2:
        raise ValueError('Unexpected original Fig.2 vector structure')
    result = {}
    for c in curves:
        points = sorted(c['pts'])
        x = 297.89883  # r_p/M=20 on right panel
        i = bisect_right([a for a, b in points], x)
        x0, y0 = points[i-1]
        x1, y1 = points[i]
        y = y0 + (y1-y0)*(x-x0)/(x1-x0)
        # Solid infinity curve is lower than dotted horizon magnitude.
        key = 'infinity' if c['top'] > 35 else 'horizon_magnitude'
        result[key] = dict(pdf_x=x, pdf_top=y,
                           plotted_flux=10**(-1-(y-27.6182)/30.36121))
    markers = json.loads(args.markers.read_text(encoding='utf-8'))
    total = sum(m['plotted_flux'] for m in markers['markers'] if m['ell'] <= 6)
    converted = total/.3**6
    report = dict(
        source='https://arxiv.org/src/2501.09806v1',
        status='cross_figure_digitization_inference_not_author_numeric_data',
        sha256_figure2=hashlib.sha256(args.figure2.read_bytes()).hexdigest(),
        sha256_figure6=markers['sha256'],
        parameters=dict(alpha=.3, rp_over_M=20, scalar_ellmax=6),
        figure2=result, figure6_sum=total,
        figure6_sum_divided_by_alpha6=converted,
        relative_difference=converted/result['infinity']['plotted_flux']-1,
        interpretation='The two graphics differ by alpha^6. Combined with epsilon^2=eta alpha^6, this supports Fig.6 per q^2 eta and Fig.2 per q^2 epsilon^2. This does not resolve the remaining solver mode discrepancy.',
        limitations='Vector interpolation and marker rounding; no formal digitization error bound. Fig.4 caption uses epsilon^2=0.1 alpha^3 whereas nearby prose specifies eta=0.1; unresolved.')
    args.output.write_bytes((json.dumps(report, indent=2)+'\n').encode('utf-8'))
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
