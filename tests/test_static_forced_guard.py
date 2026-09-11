import sys
from pathlib import Path
from types import SimpleNamespace
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from report_environment_forced_mode import argument_parser,build_metric,run


def test_cloud_resonance_is_not_sent_to_ordinary_green_solver():
    args=argument_parser().parse_args(['--metric-m','0','--scalar-ell','1'])
    cloud=SimpleNamespace(mu=.3,a=.7,rp=1.8,m=1,omega=.2)
    metric=SimpleNamespace(r0=20.,a=.7,m=0,ellmax=4,ellmin=0,omega=0.)
    with pytest.raises(ValueError,match='solvability'):
        run(args,cloud=cloud,metric=metric)


def test_static_matching_cannot_label_nonstatic_mode():
    args=argument_parser().parse_args(['--metric-m','2','--static-matching','unused.json'])
    with pytest.raises(ValueError,match='metric-m=0'):
        build_metric(args,None)
