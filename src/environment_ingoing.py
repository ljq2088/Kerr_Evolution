"""Independent horizon-regular evaluation of the threshold-cloud source.

Input h is still the covariant BL Fourier coefficient. Radial phase factors
from the advanced coordinate transformation cancel in the scalar contraction;
the cloud derivatives below retain the corresponding phase derivatives.
"""
import numpy as np
from environment_source import angular_mode


def ingoing_metric(r,theta,a):
    s=np.sin(theta)
    sigma=r*r+a*a*np.cos(theta)**2
    g=np.zeros((4,4),dtype=np.result_type(r,theta,float))
    g[0,0]=-1+2*r/sigma
    g[0,1]=g[1,0]=1
    g[0,3]=g[3,0]=-2*a*r*s*s/sigma
    g[1,3]=g[3,1]=-a*s*s
    g[2,2]=sigma
    g[3,3]=(r*r+a*a+2*a*a*r*s*s/sigma)*s*s
    return g


def ingoing_connection(r,theta,a):
    inverse=np.linalg.inv(ingoing_metric(r,theta,a))
    dg=np.zeros((4,4,4))
    dg[1]=ingoing_metric(r+1e-25j,theta,a).imag/1e-25
    dg[2]=ingoing_metric(r,theta+1e-25j,a).imag/1e-25
    gamma=np.empty((4,4,4))
    for i in range(4):
        for j in range(4):
            gamma[:,i,j]=.5*inverse@(dg[i,:,j]+dg[j,:,i]-dg[:,i,j])
    return inverse,gamma


def bl_to_ingoing_jacobian(r,a):
    delta=r*r-2*r+a*a
    jacobian=np.eye(4)
    jacobian[0,1]=-(r*r+a*a)/delta
    jacobian[3,1]=-a/delta
    return jacobian


def ingoing_cloud_hessian(cloud,r,theta):
    R,Rp=cloud.radial(r)[:,0]
    Rpp=cloud.rhs(r,[R,Rp])[1]
    S,Sp,lam=angular_mode(theta,1,1,cloud.c2)
    Spp=-np.cos(theta)/np.sin(theta)*Sp-(cloud.c2*np.cos(theta)**2-1/np.sin(theta)**2+lam)*S
    rm=2-cloud.rp
    # chi'=K_c/Delta: cancel its horizon zero analytically at threshold.
    phase_prime=cloud.omega*(r+cloud.rp)/(r-rm)
    phase_second=-2*cloud.omega/(r-rm)**2
    radial_prime=Rp+1j*phase_prime*R
    radial_second=Rpp+2j*phase_prime*Rp+(1j*phase_second-phase_prime**2)*R
    field=R*S
    grad=np.array([-1j*cloud.omega*field,radial_prime*S,R*Sp,1j*field])
    partial=np.empty((4,4),complex)
    partial[0]=-1j*cloud.omega*grad
    partial[:,0]=partial[0]
    partial[3]=1j*grad
    partial[:,3]=partial[3]
    partial[1,1]=radial_second*S
    partial[2,2]=R*Spp
    partial[1,2]=partial[2,1]=radial_prime*Sp
    inverse,gamma=ingoing_connection(r,theta,cloud.a)
    return field,partial-np.einsum('kij,k->ij',gamma,grad),inverse


def ingoing_lorenz_source(cloud,r,theta,h_bl):
    _,hessian,inverse=ingoing_cloud_hessian(cloud,r,theta)
    J=bl_to_ingoing_jacobian(r,cloud.a)
    transformed=J.T@np.asarray(h_bl)@J
    return np.einsum('ij,ij->',inverse@transformed@inverse,hessian)
