"""Validate separated source against direct tensor contraction on saved author data.
Old synchronized background is deliberately only a method test, not a Li run.
"""
import hashlib,json
from pathlib import Path
import numpy as np
from pybhpt.swsh import Yslm
from environment_source import ThresholdCloud,angular_mode
from li_separated_source import LiSeparatedSource,SPINS
from report_author_metric_source_propagation import unweight
from paper_full_tetrad import to_bl
ROOT=Path(__file__).resolve().parents[1]

def main():
    path=ROOT/'docs/environment_reproduction/sam_actual_metric_i7_h6_comparison.json'
    data=json.loads(path.read_text());cloud=ThresholdCloud();rows=[]
    assert cloud.a==data['a'] and data['m']==1
    omega=cloud.omega+1/(data['rp']**1.5+cloud.a)
    def ca(t):
        s,ds,A=angular_mode(t,1,1,cloud.c2)
        dds=-np.cos(t)/np.sin(t)*ds-(cloud.c2*np.cos(t)**2-1/np.sin(t)**2+A)*s
        return np.array([s,ds,dds])
    def ta(t):return angular_mode(t,2,2,cloud.a**2*(omega**2-cloud.mu**2))[0]
    projectors=[LiSeparatedSource(a=cloud.a,omega_c=cloud.omega,m_c=1,metric_m=1,
        spherical_lmax=3,cloud_angular=ca,target_angular=ta,quadrature=64,pmax=p) for p in (12,32)]
    x,w=np.polynomial.legendre.leggauss(64);theta=np.arccos(x)
    for row in data['rows']:
        r=row['r'];R,Rp=cloud.radial(r)[:,0];Rpp=cloud.rhs(r,[R,Rp])[1]
        for backend in ('reference','local'):
            encoded=np.asarray(row[backend]);coeff=(encoded[...,0]+1j*encoded[...,1])[:4]
            direct=0j;norm=0.
            for t,xx,ww,ss in zip(theta,x,w,ta(theta)):
                weighted=np.array([sum(coeff[j,c]*Yslm(spin,j,1,t) for j in range(max(1,abs(spin)),4)) for c,spin in enumerate(SPINS)])
                h=to_bl(unweight(weighted,r,t,cloud.a),r,t,cloud.a)
                term=2*np.pi*ww*ss*(r*r+cloud.a**2*xx*xx)*cloud.lorenz_source(r,t,h)
                direct+=term;norm+=abs(term)
            computed=[p.project(r,[R,Rp,Rpp],coeff) for p in projectors]
            enc=lambda z:[float(z.real),float(z.imag)]
            rows.append(dict(r=r,metric=backend,direct=enc(direct),separated_p12=enc(computed[0]),
                separated_p32=enc(computed[1]),p12_relative=float(abs(computed[0]/direct-1)),
                p32_relative=float(abs(computed[1]/direct-1)),angular_condition=float(norm/abs(direct))))
            print(rows[-1],flush=True)
    result=dict(status='separated_source_method_validated_on_protected_author_modes_not_Li_parameter_run',
        input=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        parameters=dict(a=cloud.a,r0=data['rp'],m_g=1,scalar_ell=2,scalar_m=2,spherical_jmax=3,q=64),
        rows=rows,limitations=['Protected spherical j<=3 only','Old synchronized cloud used solely for independent method validation',
        'This is not a=.88, full j<=18, BL-normalized Li production data'])
    (ROOT/'docs/li_alignment/separated_source_validation.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
