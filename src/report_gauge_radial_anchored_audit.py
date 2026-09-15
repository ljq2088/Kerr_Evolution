"""All 88 actual source radii: MST-anchored independent radial propagation audit."""
from pathlib import Path
import hashlib
import json
import numpy as np
from pybhpt.radial import RadialTeukolsky
from environment_mst_anchored_radial import MSTAnchoredRadial

ROOT=Path(__file__).resolve().parents[1]
A=.8771530275949366;W=1/(20**1.5+A)

def rel(x,y):return np.abs(x-y)/np.maximum(np.abs(y),1e-250)
def cx(values):return [[float(z.real),float(z.imag)] for z in values]
def packed(obj):return {bc:(obj.radialsolutions(bc),obj.radialderivatives(bc)) for bc in ('In','Up')}
def errors(values,reference,r,a):
    out={}
    for bc in ('In','Up'):
        out[bc]={}
        for j,key in enumerate(('R','Rprime')):
            e=rel(values[bc][j],reference[bc][j])
            out[bc][key]=dict(max=float(np.max(e)),radius=float(r[np.argmax(e)]),per_radius=e.tolist())
    rp=1+np.sqrt(1-a*a);delta=(r-rp)*(r-(2-rp))
    v,d=values['In'];u,ud=values['Up'];av,ad=reference['In'];au,aud=reference['Up']
    # r=20 was appended as an exact stable Wronskian anchor.
    k=int(np.flatnonzero(r==20.)[0]);power=errors.spin+1
    wr=delta**power*(v*ud-u*d);aw=delta**power*(av*aud-au*ad)
    out['Wronskian_relative_drift']=float(np.max(rel(wr,wr[k])))
    out['Wronskian_absolute_normalization_error']=float(rel(wr[k],aw[k]))
    for name,g,ag in [('horizon_kernel',u/wr[k],au/aw[k]),('infinity_kernel',v/wr[k],av/aw[k]),
            ('point_Green',np.where(r<=20,v*u[k],v[k]*u)/wr[k],np.where(r<=20,av*au[k],av[k]*au)/aw[k])]:
        e=rel(g,ag);out[name]=dict(max=float(np.max(e)),radius=float(r[np.argmax(e)]),per_radius=e.tolist())
    # Piecewise physical fields from a unit delta source; derivative jump is
    # fixed by the same W rather than separate amplitude normalization fits.
    derivative_jump=(v[k]*ud[k]-d[k]*u[k])/wr[k]
    out['unit_delta_jump_error']=float(abs(delta[k]**power*derivative_jump-1))
    return out

def main():
    cache=ROOT/'outputs/metric_cache/3e11f952b001947a3ee962e3'
    actual=np.array(sorted(float(np.load(x)['r']) for x in cache.glob('*.npz')))
    if len(actual)!=88:raise RuntimeError('Expected original actual 88-node source grid')
    r=np.unique(np.r_[actual,10.,20.]);records=[]
    variants=[(anchor,up,tol) for anchor,up in ((3.,20.),(2.1,40.)) for tol in (3e-11,3e-13,3e-14)]
    for s in (0,-1,1):
      for sign in (1,-1):
        m=sign;omega=sign*W;errors.spin=s
        native=RadialTeukolsky(s,1,m,A,omega,r);native.solve('AUTO');ref=packed(native)
        sols={};row=dict(s=s,ell=1,m=m,a=A,omega=omega,variants=[])
        for anchor,up_anchor,tol in variants:
            obj=MSTAnchoredRadial(s,1,m,A,omega,r,anchor=anchor,up_anchor=up_anchor,rtol=tol);obj.solve()
            values=packed(obj);sols[(anchor,up_anchor,tol)]=values
            comparison=errors(values,ref,r,A)
            row['variants'].append(dict(anchor=anchor,up_anchor=up_anchor,rtol=tol,diagnostics=obj.diagnostics,
                 versus_AUTO=comparison,values={bc:{'R':cx(v),'Rprime':cx(d)} for bc,(v,d) in values.items()}))
            print('MODE',s,m,'anchors',(anchor,up_anchor),'tol',tol,'In',comparison['In']['Rprime']['max'],
                'Up',comparison['Up']['Rprime']['max'],'W',comparison['Wronskian_relative_drift'],flush=True)
        row['anchor_pair_3_20_vs2p1_40_at_3e-14']=errors(sols[(3.,20.,3e-14)],sols[(2.1,40.,3e-14)],r,A)
        row['tolerance_3e-13_vs3e-14']=errors(sols[(3.,20.,3e-13)],sols[(3.,20.,3e-14)],r,A)
        row['tolerance_3e-11_vs3e-14']=errors(sols[(3.,20.,3e-11)],sols[(3.,20.,3e-14)],r,A)
        records.append(row)
    report=dict(source_grid_count=len(actual),source_grid=actual.tolist(),evaluation_grid=r.tolist(),
        source_grid_sha256=hashlib.sha256(actual.tobytes()).hexdigest(),cases=records,
        provenance={'implementation_sha256':hashlib.sha256((ROOT/'src/environment_mst_anchored_radial.py').read_bytes()).hexdigest(),
                    'report_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        scope='Massless spin0,+/-1 ell1 independent anchor/ODE reference; native radial solve only as comparator')
    path=ROOT/'docs/environment_reproduction/gauge_radial_anchored_audit.json'
    path.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(path)

if __name__=='__main__':main()
