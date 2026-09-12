import json
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from report_missing_environment_modes import build_plan


def inventory(tmp_path, modes, boundaries=None):
    data=dict(status='coverage_only_not_paper_reproduction',parameters=dict(
        alpha=.3,background='threshold-kerr',a=.877,rp=41.6,metric_ellmax=18,
        nr=8,nt=18,horizon_order=32,inner_offset=.0005,outer_source_cutoff=320,
        green_outer_by_scalar_m=boundaries or {},infinity_method_by_scalar_m={}),
        field=dict(modes=[]),infinity=dict(modes=[]),horizon=dict(modes=[
            dict(ell=l,m=m,status='missing') for l,m in modes]))
    p=tmp_path/'coverage.json';p.write_text(json.dumps(data));return p


def test_only_missing_low_dipole_uses_conjugate_branch(tmp_path):
    p=inventory(tmp_path,[(1,-1)],{'-1':4000})
    r=build_plan(p)
    assert len(r['jobs'])==1
    job=r['jobs'][0]
    assert job['metric_m']==2 and job['scalar_ells']==[] and job['conjugate_ells']==[1]
    assert '--scalar-ells' not in job['arguments']
    assert job['arguments'][job['arguments'].index('--green-outer-radius')+1]=='4000'


def test_opposite_branches_group_only_with_same_green_settings(tmp_path):
    p=inventory(tmp_path,[(3,3),(1,-1),(3,3)])
    assert len(build_plan(p)['jobs'])==1
    p=inventory(tmp_path,[(3,3),(1,-1)],{'-1':4000})
    assert len(build_plan(p)['jobs'])==2


def test_static_is_explicitly_unresolved(tmp_path):
    r=build_plan(inventory(tmp_path,[(3,1)]))
    assert r['jobs']==[] and r['unresolved_static_modes']==[[3,1]]


def test_stale_complete_inventory_does_not_hide_missing_file(tmp_path):
    p=inventory(tmp_path,[(1,-1)]);d=json.loads(p.read_text())
    d['horizon']['modes'][0].update(status='truncated_single_mode_not_converged',file='absent.json',flux={})
    p.write_text(json.dumps(d))
    assert build_plan(p)['missing_modes']==[[1,-1]]


def test_existing_response_wrong_mode_rejected(tmp_path):
    p=inventory(tmp_path,[(1,-1)]);d=json.loads(p.read_text())
    d['horizon']['modes'][0].update(status='truncated_single_mode_not_converged',file='wrong.json',flux={})
    p.write_text(json.dumps(d))
    (tmp_path/'wrong.json').write_text(json.dumps(dict(parameters=dict(scalar_ell=3,scalar_m=-1))))
    with pytest.raises(ValueError,match='mode differs'):
        build_plan(p)


def test_resolution_mismatch_not_accepted_as_completed(tmp_path):
    p=inventory(tmp_path,[(1,-1)]);d=json.loads(p.read_text())
    d['horizon']['modes'][0].update(status='truncated_single_mode_not_converged',file='mode.json',flux={})
    p.write_text(json.dumps(d))
    actual=dict(status='truncated_single_mode_not_converged',flux={},parameters=dict(
        scalar_ell=1,scalar_m=-1,alpha=.3,
        metric=dict(orbital_radius=41.6,a=.877,ellmax=18),radial_order=8,
        angular_order=10,horizon_quadrature_order=32,green_outer_radius=1000))
    response=tmp_path/'mode.json';response.write_text(json.dumps(actual))
    with pytest.raises(ValueError,match='discretization differs'):
        build_plan(p)
    actual['parameters']['angular_order']=18;response.write_text(json.dumps(actual))
    result=build_plan(p)
    assert result['missing_modes']==[] and result['verified_existing_count']==1
