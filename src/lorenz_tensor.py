"""Coordinate tensor operations for independent metric validation on Kerr."""
import numpy as np
from lorenz_corrector import tensor_covariant_derivative


def vector_covariant_derivative(g,vec):
    return [[g.partial(vec[a],b)-sum(g.gamma[c][b][a]*vec[c] for c in range(4))
             for b in range(4)] for a in range(4)]


def tensor_divergence(g,tensor):
    nabla=tensor_covariant_derivative(g,tensor)
    return [sum(g.inv[a][c]*nabla[c][a][b] for a in range(4) for c in range(4)) for b in range(4)]


def trace(g,tensor):
    return sum(g.inv[a][b]*tensor[a][b] for a in range(4) for b in range(4))


def lorenz_constraint(g,h):
    div=tensor_divergence(g,h)
    tr=trace(g,h)
    return [div[b]-g.partial(tr,b)/2 for b in range(4)]


def linearized_einstein(g,h):
    """delta Gamma -> delta Ricci -> delta Einstein, vacuum Kerr background.

    Does not assume the Lorenz condition, allowing an independent gauge test.
    """
    dh=tensor_covariant_derivative(g,h)
    dgamma=[[[sum(g.inv[c][d]*(dh[a][b][d]+dh[b][a][d]-dh[d][a][b])/2
                  for d in range(4)) for b in range(4)] for a in range(4)] for c in range(4)]
    connection_trace=[sum(dgamma[c][a][c] for c in range(4)) for a in range(4)]
    ricci=[]
    for a in range(4):
        row=[]
        for b in range(4):
            first=sum(g.partial(dgamma[c][a][b],c)
                      +sum(g.gamma[c][c][d]*dgamma[d][a][b]
                           -g.gamma[d][c][a]*dgamma[c][d][b]
                           -g.gamma[d][c][b]*dgamma[c][a][d] for d in range(4))
                      for c in range(4))
            second=g.partial(connection_trace[a],b)-sum(g.gamma[d][b][a]*connection_trace[d]
                                                        for d in range(4))
            row.append(first-second)
        ricci.append(row)
    scalar=trace(g,ricci)
    return [[ricci[a][b]-g.g[a][b]*scalar/2 for b in range(4)] for a in range(4)]


def extreme_weyl(g,h):
    """Compute delta C_lmlm and delta C_nmbnmb from coordinate curvature.

    Includes the variation of the metric lowering the Riemann index.
    Ricci subtractions vanish in these null projections; first-order tetrad
    corrections vanish because the background is type D in this tetrad.
    This does not use the reconstruction or GHP decoupling operators.
    """
    dh=tensor_covariant_derivative(g,h)
    dc=[[[sum(g.inv[c][d]*(dh[a][b][d]+dh[b][a][d]-dh[d][a][b])/2
               for d in range(4)) for b in range(4)] for a in range(4)] for c in range(4)]
    gamma=np.array([[[v.value for v in row] for row in plane] for plane in g.gamma])
    def stationary_derivative(v,index):
        return v.derivative(0 if index==1 else 1).value if index in (1,2) else 0j
    background=np.zeros((4,4,4,4),complex)
    variation=background.copy()
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    background[a,b,c,d]=stationary_derivative(g.gamma[a][b][d],c)-stationary_derivative(g.gamma[a][b][c],d)+sum(
                        gamma[a,e,c]*gamma[e,b,d]-gamma[a,e,d]*gamma[e,b,c] for e in range(4))
                    variation[a,b,c,d]=(g.partial(dc[a][b][d],c)-g.partial(dc[a][b][c],d)).value+sum(
                        gamma[a,e,c]*dc[e][b][d].value+dc[a][e][c].value*gamma[e,b,d]
                        -gamma[a,e,d]*dc[e][b][c].value-dc[a][e][d].value*gamma[e,b,c] for e in range(4))
    metric=np.array([[v.value for v in row] for row in g.g])
    perturbation=np.array([[v.value for v in row] for row in h])
    lowered=np.einsum('ae,ebcd->abcd',metric,variation)+np.einsum('ae,ebcd->abcd',perturbation,background)
    tetrad=np.array([[v.value for v in row] for row in g.tetrad])
    def project(u,v):
        return np.einsum('abcd,a,b,c,d',lowered,u,v,u,v)
    return project(tetrad[0],tetrad[2]),project(tetrad[1],tetrad[3])
