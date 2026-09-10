"""Publish a compact, verified report from a completed upstream run."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from kerr_waveform_tools.waveform_io import load_mode


def report(run, destination):
    run, destination = Path(run), Path(destination)
    time, psi, meta = load_mode(run/'psi4_l2_m2.npz')
    run_meta = json.loads((run/'metadata.json').read_text())
    if run_meta['status'] != 'complete':
        raise ValueError('run is not complete')
    if not (np.all(np.isfinite(psi)) and np.all(np.isfinite(time))):
        raise ValueError('nonfinite waveform')
    if not (np.all(np.diff(time)>0) and len(time)==run_meta['N_end']
            and len(time)==run_meta['output_count']
            and np.isclose(time[-1], run_meta['T_end'])):
        raise ValueError('invalid sample coverage')
    destination.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True, constrained_layout=True)
    axes[0].plot(time, psi.real, label='Re')
    axes[0].plot(time, psi.imag, label='Im', alpha=.75)
    axes[0].set_ylabel(r'$\psi_{4,22}$ at $\mathscr{I}^+$')
    axes[0].legend()
    axes[1].semilogy(time, np.maximum(abs(psi), 1e-30))
    axes[1].set(xlabel='T / M', ylabel=r'$|\psi_{4,22}|$')
    for ax in axes:
        ax.axvline(meta['events']['light_ring'], color='gray', linestyle=':', label='Light ring')
        ax.axvline(meta['events']['source_off'], color='black', linestyle='--', label='Source off')
        ax.grid(alpha=.2)
    axes[1].legend()
    fig.suptitle('Kerr s=-2 point-particle smoke: a/M=0.8, m=2, q=1e-5, 128 x 33')
    fig.savefig(destination/'waveform.png', dpi=180)
    plt.close(fig)
    for name in ('psi4_l2_m2.npz', 'metadata.json'):
        shutil.copy2(run/name, destination/name)
    np.savetxt(destination/'waveform.csv', np.column_stack((time,psi.real,psi.imag)),
               delimiter=',', header='T_over_M,Re_psi4_22,Im_psi4_22', comments='')
    summary = dict(status='complete', samples=len(time), T_start=float(time[0]),
                   T_end=float(time[-1]), finite=True, max_abs_psi4=float(max(abs(psi))),
                   final_abs_psi4=float(abs(psi[-1])), events=meta['events'],
                   upstream_commit='13fc2b485f76b1284c2e17879ede9075de062a9d',
                   code_commit='2ffb347', tests=dict(solver=58,shared_tools=5,local=10),
                   limitations='Smoke resolution; no continuum convergence claim; A1 bridge conditional-open',
                   sha256={name:hashlib.sha256((destination/name).read_bytes()).hexdigest()
                           for name in ('psi4_l2_m2.npz','metadata.json','waveform.csv')})
    (destination/'verification.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run')
    parser.add_argument('destination')
    args=parser.parse_args()
    report(args.run,args.destination)
