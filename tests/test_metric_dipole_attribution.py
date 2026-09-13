import json
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from report_metric_dipole_attribution import compare


def response(path,ellmax,z):
    flux=-2*abs(z)**2
    d=dict(status='truncated_single_mode_not_converged',parameters=dict(metric=dict(ellmax=ellmax),scalar_ell=0,scalar_m=0),
        samples=[dict(r=3.,weight=1.,source=[0.,0.])],z_h=[z.real,z.imag],z_inf=[0.,0.],
        flux=dict(horizon=dict(orbital_energy=flux),infinity=dict(orbital_energy=0.)))
    path.write_text(json.dumps(d));return path


def test_destructive_interference_is_retained(tmp_path):
    a=response(tmp_path/'full.json',18,.2+0j);b=response(tmp_path/'dipole.json',1,1+0j)
    row=compare(a,b)['rows']['horizon']
    assert row['dipole_flux']==pytest.approx(-2)
    assert row['higher_degree_flux']==pytest.approx(-1.28)
    assert row['interference_flux']==pytest.approx(3.2)
    assert row['full_flux']==pytest.approx(-.08)


def test_quadrature_mismatch_rejected(tmp_path):
    a=response(tmp_path/'full.json',18,1+1j);b=response(tmp_path/'dipole.json',1,1+0j)
    d=json.loads(b.read_text());d['samples'][0]['r']=4.;b.write_text(json.dumps(d))
    with pytest.raises(ValueError,match='quadrature differs'):compare(a,b)


def test_different_physics_rejected(tmp_path):
    a=response(tmp_path/'full.json',18,1+1j);b=response(tmp_path/'dipole.json',1,1+0j)
    d=json.loads(b.read_text());d['parameters']['scalar_m']=1;b.write_text(json.dumps(d))
    with pytest.raises(ValueError,match='Only metric'):compare(a,b)
