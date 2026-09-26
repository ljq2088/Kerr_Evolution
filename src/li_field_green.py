"""Real-frequency retarded Green function including negative-frequency modes."""
from types import SimpleNamespace
import numpy as np
from li_order4_green import LiOrder4Green

class LiFieldGreen(LiOrder4Green):
    def __init__(self,a,mu,omega,ell,m,**kwargs):
        if omega>=0:
            super().__init__(a,mu,omega,ell,m,**kwargs)
            return
        positive=LiOrder4Green(a,mu,-omega,ell,-m,**kwargs)
        self.__dict__.update(positive.__dict__)
        self.omega,self.m=omega,m
        self.insol=SimpleNamespace(sol=lambda r:positive.insol.sol(r).conjugate())
        self.upsol=SimpleNamespace(sol=lambda r:positive.upsol.sol(r).conjugate())
        self.w0=positive.w0.conjugate()
        self.boundary_audit=dict(positive.boundary_audit,
            conjugated_from_frequency=-omega,conjugated_from_m=-m,
            outgoing_negative_frequency=True)
