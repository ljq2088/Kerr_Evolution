"""Process-local Teukolsky method control for source-to-field comparisons.

No backend is changed on import. GSN in pybhpt 1.0.0 supports only |s|=2.
All vacuum fields and source Green functions must use the same policy.
"""
import importlib
import numpy as np
from pybhpt.radial import RadialTeukolsky as NativeRadial
from pybhpt.teuk import TeukolskyMode as NativeMode


POLICY=None
CALLS=[]


def install(spin2_method='AUTO',gauge_method='AUTO',rtol=None):
    global POLICY
    allowed={'AUTO','MST','HBL','GSN','TEUK'}
    if spin2_method not in allowed or gauge_method not in allowed-{'GSN'}:
        raise ValueError('pybhpt GSN is implemented only for spin +/-2; choose a supported gauge method')
    if rtol is not None and (not np.isfinite(rtol) or rtol<=0):raise ValueError('Positive finite rtol required')
    policy=dict(spin2_method=spin2_method,gauge_method=gauge_method,rtol=rtol)
    if POLICY is not None and POLICY!=policy:
        raise RuntimeError('Start a fresh process to change the radial policy; cached source amplitudes must not mix methods')
    if POLICY==policy:return
    # Resolve the mutually importing modules through their established entry.
    import lorenz_metric
    modules=[importlib.import_module(name) for name in ('lorenz_metric','lorenz_spin1','lorenz_spin1_chiral','lorenz_chi','lorenz_weyl')]
    for module in modules:
        for value in vars(module).values():
            if callable(getattr(value,'cache_clear',None)):value.cache_clear()
    class ControlledRadial(NativeRadial):
        diagnostic_policy=policy
        def solve(self,method='AUTO',bc=None,rtol=None):
            chosen=spin2_method if abs(self.spinweight)==2 else gauge_method
            tol=policy['rtol']
            super().solve(method=chosen,bc=bc,rtol=tol)
            for side in (('In','Up') if bc is None else (bc,)):
                values=np.r_[self.radialsolutions(side),self.radialderivatives(side)]
                if not np.all(np.isfinite(values)) or not np.any(values):
                    raise ArithmeticError(f'{chosen} returned invalid radial data: s={self.spinweight},ell={self.spheroidalmode},m={self.azimuthalmode},bc={side}')
            CALLS.append(dict(kind='radial',spin=self.spinweight,ell=self.spheroidalmode,m=self.azimuthalmode,omega=self.frequency,method=chosen,bc=bc,rtol=tol))
    class ControlledMode(NativeMode):
        diagnostic_policy=policy
        def solve(self,geo,method='AUTO',nsamples=256,teuk=None,swsh=None,tol=None):
            # Source integration tolerance differs from radial ODE tolerance;
            # preserve it. TeukolskyMode exposes method but no radial rtol.
            super().solve(geo,method=spin2_method,nsamples=nsamples,teuk=teuk,swsh=swsh,tol=tol)
            if not np.all(np.isfinite([self.amplitude('In'),self.amplitude('Up')])):
                raise ArithmeticError(f'{spin2_method} returned invalid sourced Weyl amplitudes')
            CALLS.append(dict(kind='weyl_source',method=spin2_method,radial_rtol='not_exposed_by_mode_solve'))
    for module in modules:
        if hasattr(module,'RadialTeukolsky'):module.RadialTeukolsky=ControlledRadial
    import lorenz_weyl
    lorenz_weyl.TeukolskyMode=ControlledMode
    POLICY=policy
