"""Reconstruct Re(Phi) from a saved fixed-m solution on coordinate cut planes."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors, cm
from PIL import Image
from scipy.special import sph_harm


def radial_interpolation(nodes, targets):
    weights = (-1.)**np.arange(len(nodes))
    weights[[0, -1]] *= .5
    diff = targets[:, None] - nodes[None, :]
    exact = np.isclose(diff, 0, rtol=0, atol=2e-15)
    diff[exact] = 1
    matrix = weights[None, :] / diff
    matrix /= matrix.sum(axis=1, keepdims=True)
    for row in np.flatnonzero(exact.any(axis=1)):
        matrix[row] = exact[row].astype(float)
    return matrix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', default='outputs/kerr_m2/evolution.npz')
    parser.add_argument('--output', default='outputs/kerr_m2/three_d')
    parser.add_argument('--snapshots-only', action='store_true')
    args = parser.parse_args()
    data = np.load(args.input)
    cfg = json.loads(str(data['config']))
    M, chi, m = cfg['mass'], cfg['spin'], cfg['m']
    horizon = M*(1+np.sqrt(1-chi**2))
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    r = np.linspace(horizon, 16*M, 65)
    interp = radial_interpolation(data['sigma'], horizon/r)
    # Only interpolate space; every rendered time is an existing saved time.
    modes = np.einsum('ij,tjk->tik', interp, data['u'])
    angle = np.linspace(0, 2*np.pi, 97)
    R, angle_grid = np.meshgrid(r/M, angle, indexing='ij')
    x, y = R*np.cos(angle_grid), R*np.sin(angle_grid)
    # Horizontal equatorial disk and vertical meridional disk y=0.
    theta_meridian = np.arccos(np.sin(angle))
    phi_meridian = np.where(np.cos(angle) >= 0, 0., np.pi)
    basis_eq = np.array([sph_harm(m, int(l), angle, np.pi/2) for l in data['ell']])
    basis_mer = np.array([sph_harm(m, int(l), phi_meridian, theta_meridian) for l in data['ell']])
    eq = np.einsum('trl,la->tra', modes, basis_eq).real / r[None, :, None]
    mer = np.einsum('trl,la->tra', modes, basis_mer).real / r[None, :, None]
    scale = max(float(np.max(np.abs(eq))), float(np.max(np.abs(mer))))
    norm = colors.SymLogNorm(linthresh=scale*1e-3, vmin=-scale, vmax=scale, base=10)
    cmap = plt.get_cmap('RdBu_r')
    lon, lat = np.meshgrid(np.linspace(0, 2*np.pi, 32), np.linspace(0, np.pi, 20))
    bh = horizon/M

    def draw(ax, index, labels=True):
        ax.plot_surface(x, y, np.zeros_like(x), facecolors=cmap(norm(eq[index])),
                        rstride=1, cstride=1, shade=False, linewidth=0, antialiased=False)
        ax.plot_surface(x, np.zeros_like(x), y, facecolors=cmap(norm(mer[index])),
                        rstride=1, cstride=1, shade=False, linewidth=0, antialiased=False)
        ax.plot_surface(bh*np.sin(lat)*np.cos(lon), bh*np.sin(lat)*np.sin(lon), bh*np.cos(lat),
                        color='#202631', linewidth=0, shade=True)
        ax.set(xlim=(-16,16), ylim=(-16,16), zlim=(-16,16),
               xticks=[-16,0,16], yticks=[-16,0,16], zticks=[-16,0,16])
        if labels:
            ax.set(xlabel='X / M', ylabel='Y / M', zlabel='Z / M')
        ax.set_box_aspect((1,1,1))
        ax.view_init(elev=25, azim=-55)
        ax.set_title(f'tau / M = {data["tau"][index]/M:g}', fontsize=12)
        ax.grid(alpha=.2)

    indices = [int(np.argmin(abs(data['tau']/M-t))) for t in (0,4,8,16,32,64)]
    fig = plt.figure(figsize=(15,11))
    fig.subplots_adjust(left=.02, right=.86, bottom=.08, top=.87, wspace=.12, hspace=.30)
    axes = [fig.add_subplot(2,3,i+1,projection='3d') for i in range(6)]
    for ax, idx in zip(axes,indices):
        draw(ax,idx)
    fig.suptitle(f'Kerr scalar evolution | a/M={chi}, m={m}\nRe(Phi) on equatorial and meridional coordinate cuts', fontsize=18)
    bar = fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap), cax=fig.add_axes([.9,.20,.018,.57]))
    bar.set_label('Re(Phi) | fixed symmetric-log color scale')
    fig.savefig(out/'snapshots_3d.png', dpi=160)
    plt.close(fig)
    if args.snapshots_only:
        return

    frames = []
    frame_indices = np.arange(0, len(data['tau']), 5)
    if frame_indices[-1] != len(data['tau'])-1:
        frame_indices = np.append(frame_indices, len(data['tau'])-1)
    fig = plt.figure(figsize=(8,7))
    ax = fig.add_axes([.03,.10,.79,.78], projection='3d')
    cax = fig.add_axes([.86,.22,.025,.52])
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax).set_label('Re(Phi) | fixed symlog scale')
    fig.suptitle(f'Kerr scalar field | a/M={chi}, m={m}', fontsize=17)
    fig.text(.5,.035,'Coordinate cuts: Z=0 and Y=0 | horizon <= r <= 16M\nBlack surface: r=r+ in display coordinates; not an embedding diagram',ha='center',fontsize=9)
    for count, idx in enumerate(frame_indices):
        ax.clear()
        draw(ax,int(idx))
        fig.canvas.draw()
        frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba()).copy()).convert('RGB'))
        if count % 20 == 0:
            print(f'Rendered {count+1}/{len(frame_indices)} frames',flush=True)
    frames[0].save(out/'evolution_3d.gif',save_all=True,append_images=frames[1:],duration=120,loop=0)
    plt.close(fig)
    # Verify interpolation on polynomials and endpoint exactness.
    err = float(np.max(np.abs(interp @ data['sigma']**4 - (horizon/r)**4)))
    assert err < 1e-12
    assert np.all(np.isfinite(eq)) and np.all(np.isfinite(mer))
    metadata = dict(input=str(Path(args.input).resolve()),config=cfg,
                    plotted_field='Re(sum_l u_lm Y_lm)/r; no factor of two',
                    coordinates='X=r sin(theta) cos(varphi), Y=r sin(theta) sin(varphi), Z=r cos(theta)',
                    coordinate_note='Display chart using ingoing azimuth; not Kerr-Schild Cartesian or isometric embedding',
                    radial_range_M=[horizon/M,16],cuts=['Z=0','Y=0'],
                    color_scale='fixed symmetric logarithmic',vmax=scale,linthresh=scale*1e-3,
                    frame_times_M=(data['tau'][frame_indices]/M).tolist(),
                    interpolation_polynomial_max_error=err)
    (out/'visualization.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(f'Saved 3D snapshots and animation to {out}',flush=True)


if __name__ == '__main__':
    main()
