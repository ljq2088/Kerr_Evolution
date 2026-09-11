"""Independent tensor stress integration for cloud normalization conventions.

All entries use the same stationary cloud. The r^2 measure follows the printed
paper expression in magnitude; no sign or scalar action factor is silently
changed. Only the Killing energy is used by the production source.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import simpson
from environment_source import ThresholdCloud, angular_mode


def integrate(cloud, nr, nt):
    r=cloud.rp+np.geomspace(cloud.rmin-cloud.rp,cloud.rmax-cloud.rp,nr)
    r[0],r[-1]=cloud.rmin,cloud.rmax
    x,w=np.polynomial.legendre.leggauss(nt)
    rr,xx=r[:,None],x[None,:]
    sig=rr**2+cloud.a**2*xx**2
    delta=rr**2-2*rr+cloud.a**2
    area=(rr**2+cloud.a**2)**2-cloud.a**2*delta*(1-xx**2)
    # Contravariant BL metric, evaluated independently as full arrays.
    inv=np.zeros((nr,nt,4,4))
    inv[:,:,0,0]=-area/(sig*delta)
    inv[:,:,0,3]=inv[:,:,3,0]=-2*cloud.a*rr/(sig*delta)
    inv[:,:,1,1]=delta/sig
    inv[:,:,2,2]=1/sig
    inv[:,:,3,3]=(delta-cloud.a**2*(1-xx**2))/(sig*delta*(1-xx**2))
    R,dR=cloud.radial(r)
    S,dS,_=angular_mode(np.arccos(x),1,1,cloud.c2)
    phi=R[:,None]*S
    grad=np.stack((-1j*cloud.omega*phi,dR[:,None]*S,R[:,None]*dS,1j*phi),axis=-1)
    kinetic=np.einsum('...a,...ab,...b->...',grad.conj(),inv,grad).real
    raised=np.einsum('...ab,...b->...a',inv,grad)
    # T^t_t for C=1: twice the stress tensor printed without that factor.
    minus_ttt=-(2*(raised[:,:,0]*grad[:,:,0].conj()).real-kinetic-cloud.mu**2*abs(phi)**2)
    charge_density=2*(phi.conj()*raised[:,:,0]).imag
    integrate_density=lambda f:float(2*np.pi*simpson(f@w,x=r))
    return dict(radial_points=nr,angular_points=nt,
        killing_energy=integrate_density(sig*minus_ttt),
        printed_r2_energy_magnitude_C1=integrate_density(rr**2*minus_ttt),
        noether_charge=integrate_density(sig*charge_density))


def main():
    reports=[]
    for alpha in (.2,.3):
        cloud=ThresholdCloud(alpha=alpha)
        rows=[integrate(cloud,nr,nt) for nr,nt in ((2401,32),(4801,48))]
        refined=rows[-1]
        energy=refined['killing_energy']
        q=refined['noether_charge']
        r2=refined['printed_r2_energy_magnitude_C1']
        if abs(energy-cloud.omega*q)>1e-7 or abs(energy-1)>1e-7:
            raise RuntimeError('Independent stress/charge normalization check failed')
        reports.append(dict(alpha=alpha,a=cloud.a,omega=cloud.omega,quadratures=rows,
            energy_minus_omega_charge=energy-cloud.omega*q,
            flux_multiplier_if_unit_r2_mass=energy/r2))
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/cloud_mass_measures.json'
    data=dict(status='convention_audit_not_a_source_renormalization',results=reports,
        convention='C=1 complex stress; physical Killing mass uses Sigma dr dOmega. Printed r^2 integral evaluated with positive energy sign and the same C=1.',
        limitation='Printed stress without C=1 factor would additionally divide mass by two; cannot reconcile that independently from its flux convention.')
    out.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))


if __name__=='__main__':main()
