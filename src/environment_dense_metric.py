"""Isolated deterministic angular diagnostic for complete forced-source runs.

The C++ radial solver and scalar angular projector are unchanged. Only the
four Python Lorenz angular factors use the real symmetric eigensolver.
"""
from environment_lorenz_mode import LorenzMetricMode


BACKEND = 'dense-real-evd-diagnostic'


def configure_metric_backend(dense):
    import lorenz_metric, lorenz_chi, lorenz_spin1, lorenz_spin1_chiral
    from environment_angular_diagnostic import DenseRealHarmonic, install_dense_angular_diagnostic
    modules=(lorenz_metric,lorenz_chi,lorenz_spin1,lorenz_spin1_chiral)
    installed=[module.SpinWeightedSpheroidalHarmonic is DenseRealHarmonic for module in modules]
    if any(installed) and not all(installed):
        raise RuntimeError('Mixed angular backends; start a fresh process')
    if not dense:
        if any(installed):
            raise RuntimeError('Default metric requested after dense installation; start a fresh process')
        return
    if all(installed):
        return
    caches=(lorenz_chi.chi_amplitudes,lorenz_spin1.spin1_amplitudes,
            lorenz_spin1_chiral.chiral_amplitudes,LorenzMetricMode._values)
    if any(function.cache_info().currsize for function in caches):
        raise RuntimeError('Existing metric/source amplitude caches; start a fresh process for dense diagnostics')
    install_dense_angular_diagnostic()


class DenseLorenzMetricMode(LorenzMetricMode):
    def __init__(self,*args,**kwargs):
        configure_metric_backend(True)
        super().__init__(*args,**kwargs)

    def __call__(self,r,theta):
        # Spawned workers receive a pickled object without invoking __init__.
        configure_metric_backend(True)
        return super().__call__(r,theta)

    @property
    def provenance(self):
        return dict(super().provenance,angular_backend=BACKEND)
