"""Regression and independent continued-fraction checks for the isolated nu patch."""
from pathlib import Path
import ctypes
import hashlib
import json
import math
import mpmath as mp
import numpy as np
from pybhpt.radial import RadialTeukolsky, renormalized_angular_momentum

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "outputs/pybhpt_mst_validation"
A = 0.8771530275949366
W = 1/(20**1.5+A)

def encode(z):
    z=complex(z)
    return [z.real if np.isfinite(z.real) else None, z.imag if np.isfinite(z.imag) else None]

def relative(x, y):
    return np.abs(x-y)/np.maximum(np.abs(y), 1e-290)

def peak(x):
    return float(np.max(x)) if np.isfinite(x).all() else None

def fraction(s, ell, m, a, omega, la, nu, depth):
    eps=2*mp.mpf(str(omega)); q=mp.mpf(str(a)); kap=mp.sqrt(1-q*q)
    tau=(eps-m*q)/kap; lam=mp.mpf(str(la)); j=mp.j
    def abc(n):
        v=n+nu
        al=j*eps*kap*(v+1+s+j*eps)*(v+1+s-j*eps)*(v+1+j*tau)/((v+1)*(2*v+3))
        be=-lam-s*(s+1)+v*(v+1)+eps*eps+eps*kap*tau+eps*kap*tau*(s*s+eps*eps)/(v*(v+1))
        ga=-j*eps*kap*(v-s+j*eps)*(v-s-j*eps)*(v-j*tau)/(v*(2*v-1))
        return al,be,ga
    rn=mp.mpc(0)
    for n in range(depth,0,-1):
        al,be,ga=abc(n);rn=-ga/(be+al*rn)
    ln=mp.mpc(0)
    for n in range(-depth,0):
        al,be,ga=abc(n);ln=-al/(be+ga*ln)
    al,be,ga=abc(0); terms=(be,al*rn,ga*ln)
    return sum(terms),sum(abs(x) for x in terms)

def independent_nu_check(s,ell,m,a,omega,la,nu):
    with mp.workdps(65):
        n=mp.mpc(str(nu.real),str(nu.imag))
        f=lambda v:fraction(s,ell,m,a,omega,la,v,64)[0]
        root=mp.findroot(f,(n,n+mp.mpf('1e-8')),tol=mp.mpf('1e-52'),maxsteps=30)
        r32,_=fraction(s,ell,m,a,omega,la,root,32)
        r96,_=fraction(s,ell,m,a,omega,la,root,96)
        res,scale=fraction(s,ell,m,a,omega,la,n,64)
        return {'root_65dps':str(root),'patched_nu_absolute_error':float(abs(n-root)),
                'patched_cf_absolute_residual_65dps':float(abs(res)),
                'patched_cf_relative_residual_65dps':float(abs(res)/scale),
                'root_cf_depth32_to96_residual_change':float(abs(r32-r96)),
                'root_cf_residual_depth96':float(abs(r96))}

def main():
    manifest=json.loads((BASE/'build_manifest.json').read_text())
    extension=Path(manifest['installed_extension_path'])
    assert hashlib.sha256(extension.read_bytes()).hexdigest()==manifest['installed_extension_sha256']
    lib=ctypes.CDLL(str(BASE/'libpatched_mst.so'),mode=ctypes.RTLD_LOCAL)
    f=lib.patched_mst_values
    f.argtypes=[ctypes.c_double,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_double,
                ctypes.c_double,ctypes.c_int]+[np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')]*3
    f.restype=ctypes.c_int
    cases=[(0,l,min(l,1),A,W) for l in range(7)]
    cases += [(s,l,1,A,W) for s in [-1,1] for l in [1,2,4]]
    cases += [(s,2,1,A,W) for s in [-2,2]]
    cases += [(s,1,1,0.,W) for s in [0,-1,1]]
    cases += [(s,1,-1,A,-W) for s in [0,-1,1]]
    rows=[]
    for s,ell,m,a,omega in cases:
        near=1.5 if a else 2.01
        radii=np.array([near,2.1,3.,10.,20.,40.,100.,316.],dtype=float)
        reference=RadialTeukolsky(s,ell,m,a,omega,radii);reference.solve('AUTO')
        old=renormalized_angular_momentum(s,ell,m,a,omega)
        meta=np.zeros(4);values=np.zeros((len(radii),8))
        status=f(a,s,ell,m,omega,reference.eigenvalue,len(radii),radii,meta,values)
        row={'s':s,'ell':ell,'m':m,'a':a,'omega':omega,'lambda':reference.eigenvalue,
             'status':status,'unpatched_nu':encode(old),'patched_nu':meta[:2].tolist(),
             'cpp_cf_absolute_residual':float(meta[2]),'cpp_cf_relative_residual':float(meta[3]),'radii':radii.tolist()}
        row['independent_cf']=independent_nu_check(s,ell,m,a,omega,reference.eigenvalue,complex(*meta[:2]))
        pairs=[];refpairs=[];middle=slice(1,6)
        for b,bc in enumerate(['In','Up']):
            v=values[:,b*4]+1j*values[:,b*4+1]
            d=values[:,b*4+2]+1j*values[:,b*4+3]
            av=reference.radialsolutions(bc);ad=reference.radialderivatives(bc)
            ev=relative(v,av);ed=relative(d,ad)
            row[bc]={'R':list(map(encode,v)),'R_prime':list(map(encode,d)),
                     'relative_R_error':[float(x) if np.isfinite(x) else None for x in ev],
                     'relative_R_prime_error':[float(x) if np.isfinite(x) else None for x in ed],
                     'middle_max_relative_R':peak(ev[middle]),'middle_max_relative_R_prime':peak(ed[middle]),
                     'nonfinite_radii':[float(r) for r,x,y in zip(radii,v,d) if not (np.isfinite(x) and np.isfinite(y))]}
            pairs.append((v,d));refpairs.append((av,ad))
        delta=radii*radii-2*radii+a*a
        v,d=pairs[0];u,ud=pairs[1]
        av,ad=refpairs[0];au,aud=refpairs[1]
        wr=delta**(s+1)*(v*ud-u*d);aw=delta**(s+1)*(av*aud-au*ad)
        # Green kernel connecting r=20 with the source radii; both solutions use
        # unit transmission. A delta-source measure factor is common to both.
        anchor=4;g=np.where(radii<=20,v*u[anchor],v[anchor]*u)/wr[anchor]
        ag=np.where(radii<=20,av*au[anchor],av[anchor]*au)/aw[anchor]
        row['middle_wronskian_max_relative_drift']=peak(relative(wr[middle],wr[anchor]))
        row['wronskian_at20_relative_error']=float(relative(wr[anchor],aw[anchor]))
        row['middle_green_max_relative_error']=peak(relative(g[middle],ag[middle]))
        row['middle_horizon_radiation_kernel_error']=peak(relative(u[middle]/wr[anchor],au[middle]/aw[anchor]))
        row['middle_infinity_radiation_kernel_error']=peak(relative(v[middle]/wr[anchor],av[middle]/aw[anchor]))
        row['passed_nu_regression']=(status==0 and
            row['independent_cf']['patched_nu_absolute_error']<2e-11 and
            row['independent_cf']['patched_cf_relative_residual_65dps']<1e-6)
        midok=all(row[bc]['middle_max_relative_'+q] is not None and row[bc]['middle_max_relative_'+q]<1e-7
                  for bc in ['In','Up'] for q in ['R','R_prime'])
        row['passed_bounded_regression']=(status==0 and midok and
             row['independent_cf']['patched_nu_absolute_error']<2e-11 and
             row['independent_cf']['patched_cf_relative_residual_65dps']<1e-6 and
             row['middle_green_max_relative_error'] is not None and row['middle_green_max_relative_error']<1e-7)
        rows.append(row)
        print('mode',s,ell,m,a,'pass',row['passed_bounded_regression'],
              'nuerr',row['independent_cf']['patched_nu_absolute_error'],
              'Gerror',row['middle_green_max_relative_error'],flush=True)
    report={'validation_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'build_manifest':manifest,'cases':rows,'passed':sum(x['passed_bounded_regression'] for x in rows),
            'total':len(rows),'nu_passed':sum(x['passed_nu_regression'] for x in rows),'limitations':['Only the nu translation unit is patched, and only within this diagnostic process.',
            'The remaining installed MST kernel is unchanged; direct series evaluation has separate endpoint failures.',
            'Middle-domain agreement does not establish global source-integral convergence or explain the production AUTO discrepancy.']}
    assert hashlib.sha256(extension.read_bytes()).hexdigest()==manifest['installed_extension_sha256']
    (ROOT/'docs/environment_reproduction/pybhpt_mst_patch_validation.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('passed',report['passed'],'/',report['total'])

if __name__=='__main__':main()
