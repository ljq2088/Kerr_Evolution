"""Local separated Teukolsky jets from ODE recurrence, no numerical differences."""
from lorenz_jet import Jet


def separated_jet(geometry,spin,eigenvalue,R0=1.,R1=0.,S0=1.,S1=0.):
    """Kinnersley field = R S / zeta^(|s|-s), exp(-i omega t+i m phi).

    eigenvalue is the radial lambda = A + (a omega)^2 - 2 a m omega.
    Supplying numerical In/Up and angular function values connects this local
    recurrence to actual homogeneous modes. Defaults are local test solutions.
    """
    g=geometry
    n=g.order
    r,t,a,omega,m=g.r,g.theta,g.a,g.omega,g.m
    delta=g.delta
    K=(r*r+a*a)*omega-a*m
    V=(K*K-2j*spin*(r-1)*K)/delta+4j*spin*omega*r-eigenvalue
    R=Jet(R0,n)
    R.c[1,0]=R1
    for j in range(n-1):
        rhs=-((spin+1)*2*(r-1)*R.derivative(0)+V*R)/delta
        R.c[j+2,0]=rhs.c[j,0]/((j+1)*(j+2))
    S=Jet(S0,n)
    S.c[0,1]=S1
    A=eigenvalue-a*a*omega*omega+2*a*m*omega
    U=a*a*omega*omega*t.cos()**2-2*a*omega*spin*t.cos()+spin+A-(m+spin*t.cos())**2/t.sin()**2
    for j in range(n-1):
        rhs=-t.cos()/t.sin()*S.derivative(1)-U*S
        S.c[0,j+2]=rhs.c[0,j]/((j+1)*(j+2))
    return R*S/g.zeta**(abs(spin)-spin)
