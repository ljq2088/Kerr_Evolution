"""Independent low-spin radial reference: corrected MST Cauchy data + SciPy DOP853.

No native radial solve/boundary routine is used. Unit-Teukolsky-transmission
normalization is inherited from the MST anchor and is never fitted to AUTO.
The logarithmic horizon-distance coordinate avoids cancellation in Delta.
This is a diagnostic implementation; production module bindings are unchanged.
"""
from functools import lru_cache
from pathlib import Path
import ctypes
import hashlib
import json
import numpy as np
from scipy.integrate import solve_ivp
from pybhpt.radial import RadialTeukolsky as NativeRadial

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'outputs/pybhpt_mst_validation'

@lru_cache(maxsize=1)
def _bridge():
    manifest = json.loads((BASE/'build_manifest.json').read_text())
    libpath = BASE/'libpatched_mst.so'
    if hashlib.sha256(libpath.read_bytes()).hexdigest() != manifest['library_sha256']:
        raise RuntimeError('Isolated MST bridge hash does not match its build manifest')
    lib = ctypes.CDLL(str(libpath), mode=ctypes.RTLD_LOCAL)
    f = lib.patched_mst_values
    f.argtypes = [ctypes.c_double,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_double,
        ctypes.c_double,ctypes.c_int]+[np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')]*3
    f.restype = ctypes.c_int
    return lib, f, manifest

@lru_cache(maxsize=128)
def mst_anchor(s,ell,m,a,omega,anchor):
    eigenvalue = NativeRadial(s,ell,m,a,omega,np.array([anchor])).eigenvalue
    _, f, manifest = _bridge()
    radii=np.array([anchor],dtype=float); meta=np.zeros(4); values=np.zeros((1,8))
    status=f(a,s,ell,m,omega,eigenvalue,1,radii,meta,values)
    if status or not np.isfinite(values).all():
        raise ArithmeticError(f'MST anchor invalid for {(s,ell,m,a,omega,anchor)}')
    values=values.reshape(2,4)
    data=values[:,0]+1j*values[:,1],values[:,2]+1j*values[:,3]
    return eigenvalue,data,dict(nu=meta[:2].tolist(),cf_absolute_residual=float(meta[2]),
        cf_relative_residual=float(meta[3]),anchor=float(anchor),
        mst_library_sha256=manifest['library_sha256'],normalization='unit Teukolsky transmission')

class _Propagation:
    def __init__(self,s,ell,m,a,omega,anchor,up_anchor,rtol,rmin_offset,rmax):
        self.s,self.a,self.omega,self.m=s,a,omega,m
        self.rplus=1+np.sqrt(1-a*a); self.gap=2*np.sqrt(1-a*a)
        self.eigenvalue, in_data, self.metadata=mst_anchor(s,ell,m,a,omega,anchor)
        up_eigenvalue, up_data, up_metadata=mst_anchor(s,ell,m,a,omega,up_anchor)
        if self.eigenvalue != up_eigenvalue:raise ArithmeticError('Anchor eigenvalues differ')
        v=np.array([in_data[0][0],up_data[0][1]])
        d=np.array([in_data[1][0],up_data[1][1]])
        self.anchors={'In':anchor,'Up':up_anchor}; self.rtol=rtol
        self.bounds=(self.rplus+rmin_offset,rmax)
        if not all(self.bounds[0]<point<self.bounds[1] for point in self.anchors.values()):
            raise ValueError('MST anchor must lie strictly inside propagation bounds')
        self.solutions={}; self.stats={}
        for j,bc in enumerate(('In','Up')):
            anchor=self.anchors[bc]
            x0=np.log(anchor-self.rplus)
            y0=np.array([v[j],(anchor-self.rplus)*d[j]],complex)
            scale=np.max(np.abs(y0)); y0/=scale
            halves=[]
            for side,rend in zip(('lower','upper'),self.bounds):
                solution=solve_ivp(self.rhs,(x0,np.log(rend-self.rplus)),y0,
                    method='DOP853',rtol=rtol,atol=rtol*1e-5,dense_output=True)
                if not solution.success or not np.isfinite(solution.y).all():
                    raise ArithmeticError(f'DOP853 radial propagation failed: {solution.message}')
                halves.append(solution)
                self.stats[bc+'_'+side]=dict(nfev=solution.nfev,nsteps=len(solution.t)-1)
            self.solutions[bc]=(scale,halves)
        self.metadata=dict(self.metadata,method='MST anchor + DOP853 in log(r-rplus)',
            rtol=rtol,atol=rtol*1e-5,bounds=list(self.bounds),stats=self.stats,anchor_map=self.anchors)

    def coefficients(self,r):
        q=r-self.rplus; delta=q*(q+self.gap); k=(r*r+self.a*self.a)*self.omega-self.a*self.m
        first=2*(self.s+1)*(r-1)/delta
        second=(self.eigenvalue-4j*self.s*self.omega*r-(k*k-2j*self.s*(r-1)*k)/delta)/delta
        return first,second

    def rhs(self,x,y):
        q=np.exp(x);r=self.rplus+q
        # Form Delta from q instead of subtracting nearly equal r^2-2r+a^2.
        delta=q*(q+self.gap);k=(r*r+self.a*self.a)*self.omega-self.a*self.m
        qa=2*(self.s+1)*(r-1)/(q+self.gap)
        q2b=(q*(self.eigenvalue-4j*self.s*self.omega*r)
            -(k*k-2j*self.s*(r-1)*k)/(q+self.gap))/(q+self.gap)
        return np.array([y[1],(1-qa)*y[1]+q2b*y[0]])

    def evaluate(self,bc,r):
        r=np.asarray(r,float)
        if np.any(r<self.bounds[0]) or np.any(r>self.bounds[1]):
            raise ValueError(f'Radial request outside audited propagation bounds {self.bounds}')
        if bc not in self.solutions:raise ValueError('bc must be In or Up')
        scale,halves=self.solutions[bc]
        y=np.empty((2,r.size),complex);x=np.log(r-self.rplus);mask=r<=self.anchors[bc]
        if np.any(mask):y[:,mask]=halves[0].sol(x[mask])*scale
        if np.any(~mask):y[:,~mask]=halves[1].sol(x[~mask])*scale
        v,d=y[0],y[1]/(r-self.rplus)
        first,second=self.coefficients(r)
        return v,d,-first*d+second*v

@lru_cache(maxsize=128)
def _propagation(s,ell,m,a,omega,anchor,up_anchor,rtol,rmin_offset,rmax):
    return _Propagation(s,ell,m,a,omega,anchor,up_anchor,rtol,rmin_offset,rmax)

class MSTAnchoredRadial:
    """Native RadialTeukolsky-compatible subset for diagnostic source injection.

    The independent ODE is solved once per parameter/anchor/tolerance tuple;
    repeated angular/source calls reuse dense solutions with this provenance.
    """
    default_anchor=3.0
    default_up_anchor=20.0
    default_rtol=3e-14
    default_rmin_offset=1e-5
    default_rmax=400.0

    def __init__(self,s,j,m,a,omega,r,*,anchor=None,up_anchor=None,rtol=None,rmin_offset=None,rmax=None):
        if s not in (-1,0,1) or omega==0:
            raise ValueError('This reference supports nonstatic spin 0,+/-1 modes only')
        self.spinweight=self.s=int(s);self.spheroidalmode=self.j=int(j)
        self.azimuthalmode=self.m=int(m);self.blackholespin=float(a)
        self.frequency=self.mode_frequency=self.omega=float(omega)
        self.radialpoints=np.atleast_1d(np.asarray(r,float));self.nsamples=len(self.radialpoints)
        self.anchor=self.default_anchor if anchor is None else float(anchor)
        self.up_anchor=(self.default_up_anchor if anchor is None else self.anchor) if up_anchor is None else float(up_anchor)
        self.rtol=self.default_rtol if rtol is None else float(rtol)
        self.rmin_offset=self.default_rmin_offset if rmin_offset is None else float(rmin_offset)
        self.rmax=self.default_rmax if rmax is None else float(rmax)
        self.eigenvalue=mst_anchor(s,j,m,a,omega,self.anchor)[0]
        self._values={};self.diagnostics=None

    def solve(self,method='AUTO',bc=None,rtol=None):
        if method not in ('AUTO','MST_ANCHORED'):
            raise ValueError('This diagnostic always uses MST_ANCHORED; no native fallback')
        tol=self.rtol if rtol is None else float(rtol)
        if not 2.23e-14<=tol<1e-2:raise ValueError('rtol must be between 2.23e-14 and 1e-2')
        obj=_propagation(self.s,self.j,self.m,self.blackholespin,self.omega,
            self.anchor,self.up_anchor,tol,self.rmin_offset,self.rmax)
        for side in (('In','Up') if bc is None else (bc,)):
            self._values[side]=obj.evaluate(side,self.radialpoints)
        self.diagnostics=obj.metadata
        self.provenance=dict(obj.metadata,implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())

    def _get(self,bc,deriv):
        if bc not in self._values:raise RuntimeError('Call solve before retrieving radial data')
        return self._values[bc][deriv].copy()
    def radialpoint(self,pos):return self.radialpoints[pos]
    def radialsolutions(self,bc):return self._get(bc,0)
    def radialderivatives(self,bc):return self._get(bc,1)
    def radialderivatives2(self,bc):return self._get(bc,2)
    def radialsolution(self,bc,pos):return self._get(bc,0)[pos]
    def radialderivative(self,bc,pos):return self._get(bc,1)[pos]
    def radialderivative2(self,bc,pos):return self._get(bc,2)[pos]
    def __call__(self,bc,deriv=0):return self._get(bc,deriv)
