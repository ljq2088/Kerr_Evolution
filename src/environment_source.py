"""Normalized massive cloud and covariant Lorenz-source contraction in BL.

Coordinates (t,r,theta,phi), M=1, signature -+++, complex scalar C=1.
This module does not invent a point-particle metric: callers must supply h.
"""
import numpy as np
from scipy.integrate import solve_ivp, simpson
from scipy.special import lpmv, gammaln
from environment_cloud import cloud_211, angular_eigenvalue, radial_coefficients


def kerr_metric(r, theta, a):
    s, c = np.sin(theta), np.cos(theta)
    sigma, delta = r*r+a*a*c*c, r*r-2*r+a*a
    g = np.zeros((4, 4), dtype=np.result_type(r, theta, float))
    g[0, 0] = -1+2*r/sigma
    g[0, 3] = g[3, 0] = -2*a*r*s*s/sigma
    g[1, 1], g[2, 2] = sigma/delta, sigma
    g[3, 3] = (r*r+a*a+2*a*a*r*s*s/sigma)*s*s
    return g


def connection(r, theta, a):
    """Complex-step derivatives of the analytic BL metric (off the axes)."""
    inverse = np.linalg.inv(kerr_metric(r, theta, a))
    dg = np.zeros((4, 4, 4))  # derivative index first
    dg[1] = kerr_metric(r+1e-25j, theta, a).imag/1e-25
    dg[2] = kerr_metric(r, theta+1e-25j, a).imag/1e-25
    gamma = np.empty((4, 4, 4))
    for i in range(4):
        for j in range(4):
            gamma[:, i, j] = .5*inverse@(dg[i, :, j]+dg[j, :, i]-dg[:, i, j])
    return inverse, gamma


def angular_mode(theta, ell, m, c2, size=20):
    """S and dS/dtheta, with integral |S exp(im phi)|² dOmega = 1."""
    if m < 0:
        s, ds, lam = angular_mode(theta, ell, -m, c2, size)
        return (-1)**abs(m)*s, (-1)**abs(m)*ds, lam
    if ell < m:
        raise ValueError('Require ell >= |m|')
    ls = np.arange(m, m+size+1)
    cosine = np.zeros((size+1, size+1))
    for j, l in enumerate(ls[:-1]):
        cosine[j, j+1] = cosine[j+1, j] = np.sqrt(((l+1)**2-m*m)/((2*l+1)*(2*l+3)))
    matrix = np.diag(ls[:size]*(ls[:size]+1.))-c2*(cosine@cosine)[:size, :size]
    eigenvalues, vectors = np.linalg.eigh(matrix)
    coeff = vectors[:, ell-m]
    coeff *= np.sign(coeff[ell-m])
    theta = np.asarray(theta)
    x, st = np.cos(theta), np.sin(theta)
    value, derivative = np.zeros_like(x), np.zeros_like(x)
    for l, b in zip(ls[:size], coeff):
        norm = np.sqrt((2*l+1)/(4*np.pi)*np.exp(gammaln(l-m+1)-gammaln(l+m+1)))
        p = lpmv(m, l, x)
        pm = lpmv(m, l-1, x) if l > m else np.zeros_like(x)
        value += b*norm*p
        derivative += b*norm*(l*x*p-(l+m)*pm)/st
    return value, derivative, eigenvalues[ell-m]


class ThresholdCloud:
    """|211>, physical field normalized to requested Killing mass."""
    def __init__(self, mass=1., radial_points=2401, alpha=.3):
        result = cloud_211(outer_efolds=45., horizon_offset=1e-6, rtol=2e-11,mu=alpha)
        self.a, self.omega, self.rp = (result[k] for k in
                                      ('a_over_M', 'M_omega_c', 'r_plus_over_M'))
        self.mu, self.m = alpha, 1
        self.c2 = self.a**2*(self.omega**2-self.mu**2)
        self.lam = angular_eigenvalue(1, 1, self.c2)
        self.rmin = self.rp+1e-6
        k = np.sqrt(self.mu**2-self.omega**2)
        self.rmax, self.match = 45/k, 2/self.mu**2
        def rhs(r, state):
            d, dp, v = radial_coefficients(r, self.a, self.mu, self.omega, 1, self.lam)
            return [state[1], -(dp*state[1]+v*state[0])/d]
        self.rhs = rhs
        vh = -self.mu**2*self.rp**2-self.a**2*self.omega**2+2*self.a*self.omega-self.lam
        slope = -vh/(2*(self.rp-1))
        self.left = solve_ivp(rhs, (self.rmin, self.match), [1+slope*1e-6, slope],
                              dense_output=True, rtol=2e-11, atol=2e-13)
        beta = (2*self.omega**2-self.mu**2)/k-1
        self.right = solve_ivp(rhs, (self.rmax, self.match), [1., -k+beta/self.rmax],
                               dense_output=True, rtol=2e-11, atol=2e-13)
        if not self.left.success or not self.right.success:
            raise RuntimeError('Cloud profile integration failed')
        self.glue = self.left.y[0, -1]/self.right.y[0, -1]
        self.amplitude = 1.
        energy, charge = self.integrals(radial_points)
        if mass <= 0 or energy <= 0:
            raise ValueError('Mass and cloud energy must be positive')
        self.amplitude = np.sqrt(mass/energy)
        self.mass, self.charge = mass, charge*self.amplitude**2

    def radial(self, r):
        r = np.atleast_1d(r).astype(float)
        if np.any(r < self.rmin) or np.any(r > self.rmax):
            raise ValueError('Requested radius outside integrated cloud domain')
        out = np.empty((2, r.size))
        mask = r <= self.match
        if mask.any():
            out[:, mask] = self.left.sol(r[mask])
        if (~mask).any():
            out[:, ~mask] = self.glue*self.right.sol(r[~mask])
        return self.amplitude*out

    def integrals(self, radial_points=2401):
        """Direct stress-energy and Noether charge integrals, not E=omega Q."""
        r = self.rp+np.geomspace(self.rmin-self.rp, self.rmax-self.rp, radial_points)
        r[0], r[-1] = self.rmin, self.rmax
        R, dR = self.radial(r)
        x, weights = np.polynomial.legendre.leggauss(40)
        s, ds, _ = angular_mode(np.arccos(x), 1, 1, self.c2)
        rr, xx = r[:, None], x[None, :]
        sigma, delta = rr**2+self.a**2*xx**2, rr**2-2*rr+self.a**2
        gtt = -((rr**2+self.a**2)**2-self.a**2*delta*(1-xx**2))/(sigma*delta)
        gtp = -2*self.a*rr/(sigma*delta)
        gpp = (delta-self.a**2*(1-xx**2))/(sigma*delta*(1-xx**2))
        field2 = R[:, None]**2*s**2
        e = sigma*((-gtt*self.omega**2+gpp+self.mu**2)*field2
                   +delta/sigma*dR[:, None]**2*s**2+R[:, None]**2*ds**2/sigma)
        q = 2*sigma*(-gtt*self.omega+gtp)*field2
        return (2*np.pi*simpson(e@weights, x=r), 2*np.pi*simpson(q@weights, x=r))

    def hessian(self, r, theta):
        """phi0 and covariant Hessian at t=phi=0; same harmonic factors omitted."""
        R, dR = self.radial(r)[:, 0]
        ddR = self.rhs(r, [R, dR])[1]
        S, dS, lam = angular_mode(np.asarray(theta), 1, 1, self.c2)
        ddS = -np.cos(theta)/np.sin(theta)*dS-(self.c2*np.cos(theta)**2
                                                      -1/np.sin(theta)**2+lam)*S
        field = R*S
        grad = np.array([-1j*self.omega*field, dR*S, R*dS, 1j*field])
        partial = np.empty((4, 4), complex)
        partial[0] = -1j*self.omega*grad
        partial[:, 0] = partial[0]
        partial[3] = 1j*grad
        partial[:, 3] = partial[3]
        partial[1, 1], partial[2, 2] = ddR*S, R*ddS
        partial[1, 2] = partial[2, 1] = dR*dS
        inverse, gamma = connection(r, theta, self.a)
        return field, partial-np.einsum('kij,k->ij', gamma, grad), inverse

    def lorenz_source(self, r, theta, h_covariant):
        """h^{ab} Hessian_ab. Caller must establish Lorenz gauge and provenance."""
        _, hessian, inverse = self.hessian(r, theta)
        h = np.asarray(h_covariant, dtype=complex)
        if h.shape != (4, 4) or not np.allclose(h, h.T):
            raise ValueError('Expected symmetric covariant BL metric perturbation')
        return np.einsum('ij,ij->', inverse@h@inverse, hessian)


def project_source(cloud, radii, orbital_radius, ell, m, metric_mode, ntheta=48):
    """Project Sigma h:H onto S_lm; metric_mode(r,theta) is the m_g=m-1 mode.

    Input must be covariant BL Lorenz h, with q stripped. The callback returns
    its complex Fourier coefficient, not its real part or trace reverse.
    No factor of 2 is added: a real metric's negative-m_g coefficient is its
    complex conjugate. The cloud mass convention propagates into J unchanged.
    """
    omega_p = 1/(orbital_radius**1.5+cloud.a)
    omega = cloud.omega+(m-cloud.m)*omega_p
    x, weights = np.polynomial.legendre.leggauss(ntheta)
    theta = np.arccos(x)
    s = angular_mode(theta, ell, m, cloud.a**2*(omega**2-cloud.mu**2))[0]
    projected = []
    for r in radii:
        source = np.array([cloud.lorenz_source(r, t, metric_mode(r, t)) for t in theta])
        projected.append(2*np.pi*np.dot(weights, s*(r*r+cloud.a**2*x*x)*source))
    return omega, np.asarray(projected)
