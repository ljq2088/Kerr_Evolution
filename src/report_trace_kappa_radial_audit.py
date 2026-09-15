"""Actual dipole s=0 trace/kappa radial and bounded source-node audit.

No production module is edited. This extends the prior kappa-only full-grid
comparison by changing trace and kappa together, testing finite difference
step/solver controls, and comparing trace against independent pybhpt backends.
"""
import json,hashlib,time
from pathlib import Path
import numpy as np
from pybhpt.radial import RadialTeukolsky
import lorenz_kappa as original
import lorenz_metric as metric_module
import environment_radial_variation as rv
import environment_trace_variation as tv
from lorenz_ghp import KerrGHP
from lorenz_spin1 import cky_tensor
from lorenz_tensor import vector_covariant_derivative
from environment_source import ThresholdCloud,angular_mode,kerr_metric

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'docs/environment_reproduction/trace_kappa_radial_precision_audit.json'
A=.8771530275949366;R0=20.;ELL=M=1;OMEGA=1/(R0**1.5+A)

def enc(x):
    x=np.asarray(x,complex)
    return np.stack([x.real,x.imag],axis=-1).tolist()

def refined_variation():
    old_class=tv.RadialMassVariation;old_h=rv.horizon_mass_boundary;old_i=rv.infinity_mass_boundary
    def cls(*args,**kwargs):return old_class(*args,offset=1e-5,**kwargs)
    def horizon(*args,**kwargs):kwargs['order']=8;return old_h(*args,**kwargs)
    def infinity(*args,**kwargs):kwargs['order']=12;return old_i(*args,**kwargs)
    tv.RadialMassVariation=cls;rv.horizon_mass_boundary=horizon;rv.infinity_mass_boundary=infinity
    try:return tv.TraceMassVariation(R0,A,ELL,M,rmax=8000.,rtol=5e-14)
    finally:tv.RadialMassVariation=old_class;rv.horizon_mass_boundary=old_h;rv.infinity_mass_boundary=old_i

def correction_metric(g,trace,kappa):
    ky=cky_tensor(g)
    xi=[sum(ky[i][j]*g.inv[j][k]*g.partial(trace,k)/2 for j in range(4) for k in range(4))
         +g.partial(kappa,i) for i in range(4)]
    d=vector_covariant_derivative(g,xi)
    return np.array([[-(d[i][j].value+d[j][i].value)*1j/OMEGA for j in range(4)] for i in range(4)])

def jet_components(x):
    indices=[(0,0),(1,0),(0,1),(2,0),(1,1),(0,2),(3,0),(4,0)]
    return np.array([x.derivative_value(*ij) for ij in indices])

def main():
    file=ROOT/'docs/environment_reproduction/forced_mode_nr8_nt18_L1_mg-1_sl0_inner0.0005_outer320_log_h32.json'
    raw=file.read_bytes();data=json.loads(raw)
    indices=[0,16,31,35,43,55,63,79,87]
    samples=[data['samples'][i] for i in indices];radii=np.array([s['r'] for s in samples])
    # Build new diagnostic variations without changing source-producing defaults.
    default_var=tv.TraceMassVariation(R0,A,ELL,M,rmax=2000.,rtol=1e-12)
    refined=refined_variation()
    step=.01*OMEGA**2
    def pair(g,name):
        if name=='baseline':return original.trace_field_jet(g,R0,ELL),original.kappa_jet(g,R0,ELL)
        if name=='half_step':return original.trace_field_jet(g,R0,ELL),original.kappa_jet(g,R0,ELL,step=step/2)
        if name=='finite_refined':return original.trace_field_jet(g,R0,ELL,rtol=5e-14,rmax=4000.),original.kappa_jet(g,R0,ELL,rtol=5e-14,rmax=4000.)
        var=default_var if name=='analytic' else refined
        h,dh=var.field_jet(g);return h,-1j*OMEGA*dh
    names=['baseline','half_step','finite_refined','analytic','analytic_refined']
    result=dict(status='running_bounded_trace_kappa_radial_audit',parameters=dict(a=A,r0=R0,ell=ELL,m=M,omega=OMEGA,step=step,angular_order=18,source_indices=indices),
       input_file=file.name,input_sha256=hashlib.sha256(raw).hexdigest(),rows=[],trace_backends=[],
       configurations=dict(baseline='Original central Richardson difference, order4/6 endpoints, rmax2000, offset1e-4, rtol1e-12',
          half_step='Original settings, halve auxiliary mass-squared difference step',
          finite_refined='Finite derivative with rtol5e-14 and rmax4000',
          analytic='Differentiated radial/angular ODE at original cutoffs',
          analytic_refined='Differentiated ODE: rmax8000, offset1e-5, horizon order8, infinity order12, rtol5e-14'),
       source_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [Path(__file__),ROOT/'src/lorenz_kappa.py',ROOT/'src/lorenz_mode_jet.py',ROOT/'src/environment_trace_variation.py',ROOT/'src/environment_radial_variation.py']},
       limitations=['Nine actual radial nodes, angular quadrature order18; not a new complete radial integral.',
       'Other metric sectors, compact chi and cloud remain fixed; source arrays are historical and hashed.',
       'Independent pybhpt scalar Green checks share the physical source and separation conventions.',
       'Existing full-grid kappa-only evidence is cited separately, not relabeled as new.'])
    def save():
        tmp=OUT.with_suffix('.tmp');tmp.write_text(json.dumps(result,indent=2)+'\n');tmp.replace(OUT)
    source_eq=-16*np.pi*angular_mode(np.pi/2,ELL,M,A*A*OMEGA*OMEGA)[0]
    metric=kerr_metric(R0,np.pi/2,A);ut=1/np.sqrt(-metric[0,0]-2*OMEGA*metric[0,3]-OMEGA**2*metric[3,3]);source_eq/=ut
    query=np.sort(np.concatenate([radii,[R0]]));j0=int(np.where(query==R0)[0][0])
    radial,zi,zh=original.scalar_resolvent(R0,A,ELL,M,0.)
    reference=np.array([(zh*radial.insol.sol(r) if r<R0 else zi*radial.upsol.sol(r)) for r in query])
    for method in ['AUTO','HBL','TEUK','MST']:
        t=time.perf_counter();mode=RadialTeukolsky(0,ELL,M,A,OMEGA,query);mode.solve(method=method)
        u=np.array([[mode.radialsolution('In',j),mode.radialderivative('In',j)] for j in range(len(query))])
        v=np.array([[mode.radialsolution('Up',j),mode.radialderivative('Up',j)] for j in range(len(query))])
        if not np.isfinite(u).all() or not np.isfinite(v).all():
            row=dict(method=method,status='backend_nonfinite_radial_output',
                nonfinite_In_entries=int((~np.isfinite(u)).sum()),nonfinite_Up_entries=int((~np.isfinite(v)).sum()),
                lambda_radial=mode.eigenvalue,seconds=time.perf_counter()-t)
            result['trace_backends'].append(row);save();print('trace backend',method,row['status'],flush=True);continue
        W=(R0*R0-2*R0+A*A)*(u[j0,0]*v[j0,1]-v[j0,0]*u[j0,1])
        actual=np.array([(u[j]*v[j0,0] if r<R0 else u[j0,0]*v[j])*source_eq/W for j,r in enumerate(query)])
        row=dict(method=method,status='finite_scalar_trace_comparison',lambda_radial=mode.eigenvalue,lambda_angular=radial.lam,
           eigenvalue_convention_error=float(abs(mode.eigenvalue-(radial.lam+A*A*OMEGA*OMEGA-2*A*M*OMEGA))),
           radial_trace=enc(actual),reference_trace=enc(reference),radii=query.tolist(),
           max_relative_value_error=float(max(abs(actual[:,0]/reference[:,0]-1))),
           max_relative_derivative_error=float(max(abs(actual[:,1]/reference[:,1]-1))),seconds=time.perf_counter()-t)
        result['trace_backends'].append(row);save();print('trace backend',method,row['max_relative_value_error'],flush=True)
    cloud=ThresholdCloud(alpha=.3);assert abs(cloud.a-A)<1e-13
    x,w=np.polynomial.legendre.leggauss(18);theta=np.arccos(x)
    spherical=angular_mode(theta,0,0,A*A*(data['parameters']['omega']**2-.3**2))[0]
    for index,r in zip(indices,radii):
        t0=time.perf_counter();g=KerrGHP(float(r),1.1,A,omega=OMEGA,m=1,order=4)
        base=pair(g,'baseline');local={}
        for name in names:
            jets=pair(g,name);hh,kk=map(jet_components,jets);bh,bk=map(jet_components,base)
            local[name]=dict(trace_derivatives=enc(hh),kappa_derivatives=enc(kk),
                max_relative_trace_derivatives=float(max(abs(hh/bh-1))),max_relative_kappa_derivatives=float(max(abs(kk/bk-1))))
        delta={name:0j for name in names if name!='baseline'}
        fresh_source=0j
        for tt,xx,ww,ss in zip(theta,x,w,spherical):
            gg=KerrGHP(float(r),float(tt),A,omega=OMEGA,m=1,order=2)
            hb,kb=pair(gg,'baseline');_,hess,inv=cloud.hessian(r,tt);raised=inv@hess@inv
            factor=2*np.pi*ww*ss*(r*r+A*A*xx*xx)
            metric_jets=metric_module.nonstatic_metric(float(r),float(tt),R0,A,ELL,M,order=6)[1]
            metric=np.array([[v.value for v in row] for row in metric_jets],complex).conjugate()
            fresh_source+=factor*np.einsum('ij,ij->',metric,raised)
            for name in delta:
                hh,kk=pair(gg,name);dh=correction_metric(gg,hh-hb,kk-kb).conjugate()
                delta[name]+=factor*np.einsum('ij,ij->',dh,raised)
        stored_source=complex(*data['samples'][index]['source'])
        base_source=fresh_source
        row=dict(index=index,r=r,local_at_theta1p1=local,baseline_source=enc(base_source),stored_source=enc(stored_source),
            fresh_vs_stored_relative_change=float(abs(base_source/stored_source-1)),
            projected_source={name:dict(delta=enc(dz),new=enc(base_source+dz),relative_complex_change=float(abs(dz/base_source))) for name,dz in delta.items()},seconds=time.perf_counter()-t0)
        result['rows'].append(row);save();print('source node',index,r,{n:v['relative_complex_change'] for n,v in row['projected_source'].items()},flush=True)
    result['status']='completed_bounded_trace_kappa_radial_audit';save()

if __name__=='__main__':main()
