import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from environment_source import ThresholdCloud,kerr_metric
from environment_ingoing import ingoing_metric,ingoing_cloud_hessian,bl_to_ingoing_jacobian,ingoing_lorenz_source


def test_ingoing_metric_hessian_and_source_transform():
    cloud=ThresholdCloud(alpha=.3)
    rng=np.random.default_rng(614)
    tensor=rng.normal(size=(4,4))
    tensor=(tensor+tensor.T)/2
    for r in (3.,10.,30.):
        J=bl_to_ingoing_jacobian(r,cloud.a)
        np.testing.assert_allclose(J.T@kerr_metric(r,1.1,cloud.a)@J,ingoing_metric(r,1.1,cloud.a),atol=1e-12)
        _,H,_=cloud.hessian(r,1.1)
        _,regular,_=ingoing_cloud_hessian(cloud,r,1.1)
        np.testing.assert_allclose(J.T@H@J,regular,rtol=1e-10,atol=1e-13)
        np.testing.assert_allclose(ingoing_lorenz_source(cloud,r,1.1,tensor),cloud.lorenz_source(r,1.1,tensor),rtol=1e-10,atol=1e-13)
    for offset in (1e-2,1e-3,1e-4):
        field,H,inverse=ingoing_cloud_hessian(cloud,cloud.rp+offset,1.1)
        np.testing.assert_allclose(np.einsum('ij,ij->',inverse,H),cloud.mu**2*field,rtol=1e-9,atol=1e-12)
