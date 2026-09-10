"""Generate independently computed background artifacts; no forced wake claimed."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from environment_source import ThresholdCloud, angular_mode


def main(alpha=.3):
    out = Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    if alpha!=.3:
        out=out/'alpha_02'
    out.mkdir(exist_ok=True,parents=True)
    cloud = ThresholdCloud(alpha=alpha)
    r = np.geomspace(cloud.rmin, cloud.rmax, 1601)
    R, dR = cloud.radial(r)
    e, q = cloud.integrals(4801)
    residual = []
    for rr in (cloud.rp+.01,3.5,20.,41.6,100.):
        for t in (.3,1.2,2.6):
            f,h,inv = cloud.hessian(rr,t)
            residual.append(abs(np.einsum('ij,ij',inv,h)/(cloud.mu**2*f)-1))
    report = dict(status='background_only_not_forced_environment', a=cloud.a,
                  mu=cloud.mu, omega_c=cloud.omega, normalized_killing_mass=1.,
                  direct_energy=float(e), direct_charge=float(q),
                  energy_minus_omega_charge=float(e-cloud.omega*q),
                  maximum_relative_KG_trace_residual=float(max(residual)),
                  field_convention='C=1 complex scalar; integral |S exp(im phi)|² dOmega=1',
                  paper_phi10_scale='multiply this mass=1 field by alpha**(-3)')
    (out/'background_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez_compressed(out/'background_211.npz', r=r, R=R, dR=dR,
                        a=cloud.a, omega=cloud.omega, mu=cloud.mu, mass=1.)
    fig, axes = plt.subplots(1,2,figsize=(10,4), constrained_layout=True)
    axes[0].plot(r,R)
    axes[0].set(xlim=(cloud.rp,150),xlabel='Boyer-Lindquist r/M',ylabel='R(r), cloud mass = 1',
                title=f'Kerr threshold cloud |211>, alpha = {alpha}')
    x=np.linspace(-80,80,301)
    xx,yy=np.meshgrid(x,x)
    rr=np.hypot(xx,yy)
    vals=np.full_like(rr,np.nan)
    mask=rr>=cloud.rmin
    vals[mask]=abs(cloud.radial(rr[mask])[0]*angular_mode(np.pi/2,1,1,cloud.c2)[0])
    im=axes[1].pcolormesh(xx,yy,vals,shading='auto',cmap='magma')
    axes[1].set(aspect='equal',xlabel='r cos(phi)/M',ylabel='r sin(phi)/M',
                title='Background |Phi|, equatorial coordinate slice')
    fig.colorbar(im,ax=axes[1],label='|Phi| (mass = 1)')
    fig.savefig(out/'background_211.png',dpi=160)
    plt.close(fig)
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--alpha',type=float,choices=(.2,.3),default=.3)
    main(parser.parse_args().alpha)
