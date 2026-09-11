"""Record completed high-L channels without treating missing channels as zero."""
import hashlib
import json
from pathlib import Path


def main():
    directory = Path('docs/environment_reproduction')
    inputs = {}

    def read(name):
        raw = (directory / name).read_bytes()
        inputs[name] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)

    markers = read('paper_v1_flux_markers.json')
    reference = {(row['ell'], row['m']): row['plotted_flux']
                 for row in markers['markers']}
    comparison = []
    for name in ('scalar_batch_mg1_e732460f3ebb.json',
                 'scalar_batch_mg2_55f1591fe654.json'):
        batch = read(name)
        if batch['status'] != 'batch_completed_finite_resolution_not_converged':
            raise ValueError(f'Incomplete batch: {name}')
        for channel in batch['channels']:
            key = channel['scalar_ell'], channel['scalar_m']
            if key not in reference:
                continue
            flux = channel['flux']['infinity']['orbital_energy']
            comparison.append(dict(ell=key[0], m=key[1], source=channel['file'],
                computed=flux, digitized=reference[key], ratio=flux/reference[key]))

    stem = 'forced_mode_nr8_nt18_L18_mg1_sl12_gh0.0001_go'
    suffix = '_inner0.0005_outer320_log_h32.json'
    near = read(stem+'4000'+suffix)
    far = read(stem+'8000'+suffix)
    if [(x['r'], x['source']) for x in near['samples']] != [
            (x['r'], x['source']) for x in far['samples']]:
        raise ValueError('Boundary audit must reuse identical sources')
    boundary = dict(scalar_ell=12, scalar_m=2, outer_radii=[4000, 8000],
        reused_samples=sum(x['reused'] for x in far['samples']),
        relative_flux_changes={boundary:
            far['flux'][boundary]['orbital_energy']/near['flux'][boundary]['orbital_energy']-1
            for boundary in ('infinity', 'horizon')})
    result = dict(status='partial_high_L_comparison_not_full_figure6',
        inputs_sha256=inputs, comparisons=comparison, boundary_audit=boundary,
        limitations=[
            'Only completed m=2 and m=3 infinity sectors compared; no missing-mode sum',
            'Paper values digitized from vector figure, not author numerical data',
            'Units inferred as q^2 Mc/M by independent Figures 2 and 6 closure',
            'Outer-boundary audit does not establish source-grid or metric convergence',
            'High-L metric matching accuracy remains separately unresolved'])
    output = directory/'rp20_high_L_partial_audit.json'
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(output)


if __name__ == '__main__':
    main()
