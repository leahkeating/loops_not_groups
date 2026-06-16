# negative binomial pgfs

import numpy as np

# Negative-binomial networks: the single-edge count and the triangle count are
#   each negative-binomial distributed with dispersion r (here r = 3). A node's
#   total mean degree is fixed at `deg`, split as `mean_edges = deg - 2*mean_tri`
#   single-edge stubs and `mean_tri` triangles. Writing p = r/(r+mean), each NB
#   PGF is (p/(1-(1-p)z))**r, so G0 is the product of the edge and triangle PGFs.
#   As r -> inf this recovers the Poisson case in ../poisson_pgfs.py.
def G0(x, y, mean_tri, deg=6):
    r = 3
    mean_edges = deg - 2*mean_tri
    p_s = r/(r+mean_edges)
    p_t = r/(r+mean_tri)
    return (((p_s*p_t)/((1-(1-p_s)*x)*(1-(1-p_t)*y)))**r)

def Gs(x, y, mean_tri, deg=6):
    r = 3
    mean_edges = deg - 2*mean_tri
    p_s = r/(r+mean_edges)
    return p_s*G0(x, y, mean_tri, deg)/(1-(1-p_s)*x)

def Gt(x, y, mean_tri, deg=6):
    r = 3
    p_t = r/(r+mean_tri)
    return p_t*G0(x, y, mean_tri, deg)/(1-(1-p_t)*y)
