"""Backend isolation must hold in fresh processes and spawned cache workers."""
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]


def test_dense_source_backend_spawn_and_cache_isolation(tmp_path):
    script=tmp_path/'check_dense.py'
    script.write_text('''
from pathlib import Path
import sys
import numpy as np
from environment_lorenz_mode import LorenzMetricMode
from environment_dense_metric import DenseLorenzMetricMode,configure_metric_backend,BACKEND
from environment_metric_sampling import precompute_metric

def main():
    folder=Path(sys.argv[1])
    base=LorenzMetricMode(6.,.6,2,2)
    original,audit=precompute_metric(base,[4.5],[.7],folder,workers=1)
    # Original sampling ran in a child; the parent has no source caches.
    dense=DenseLorenzMetricMode(6.,.6,2,2)
    sampled,new=precompute_metric(dense,[4.5],[.7],folder,workers=1)
    assert new['cache_key']!=audit['cache_key']
    assert new['reconstructed_radii']==1
    assert sampled.provenance['angular_backend']==BACKEND
    np.testing.assert_allclose(sampled(4.5,.7),dense(4.5,.7),rtol=1e-12,atol=1e-14)
    again,resumed=precompute_metric(dense,[4.5],[.7],folder,workers=1)
    assert resumed['cached_radii']==1 and resumed['reconstructed_radii']==0
    np.testing.assert_array_equal(again(4.5,.7),sampled(4.5,.7))
    try:
        configure_metric_backend(False)
    except RuntimeError:
        pass
    else:
        raise AssertionError('Default backend accepted after dense installation')

if __name__=='__main__':main()
''')
    environment=dict(os.environ,PYTHONPATH=str(ROOT/'src'),OPENBLAS_NUM_THREADS='1')
    subprocess.run([sys.executable,str(script),str(tmp_path/'cache')],env=environment,check=True,timeout=180)


def test_dense_rejects_existing_default_amplitudes():
    code='''
from lorenz_chi import chi_amplitudes
from environment_dense_metric import configure_metric_backend
chi_amplitudes(6.,.6,2,2)
try:
    configure_metric_backend(True)
except RuntimeError as error:
    assert 'caches' in str(error)
else:
    raise AssertionError('Reused incompatible source cache')
'''
    environment=dict(os.environ,PYTHONPATH=str(ROOT/'src'),OPENBLAS_NUM_THREADS='1')
    subprocess.run([sys.executable,'-c',code],env=environment,check=True,timeout=180)
