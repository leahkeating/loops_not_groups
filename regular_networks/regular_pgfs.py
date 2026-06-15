# regular pgfs

import numpy as np

# Regular networks: every node has the same fixed degree `deg` and the same
#   number of triangles `num_tri`. Each triangle uses up 2 of the node's stubs,
#   leaving (deg - 2*num_tri) single edges. Because the distribution p_st is
#   deterministic, the generating functions are single monomials rather than
#   the Poisson sums in poisson_pgfs.py.
def G0(x, y, num_tri, deg=6):
    x_exponent = np.max(deg - 2*num_tri, 0)
    y_exponent = np.max(num_tri, 0)
    return x**x_exponent * y**y_exponent

def Gs(x, y, num_tri, deg=6):
    x_exponent = np.max(deg - 2*num_tri - 1, 0)
    y_exponent = np.max(num_tri, 0)
    return x**x_exponent * y**y_exponent

def Gt(x, y, num_tri, deg=6):
    x_exponent = np.max(deg - 2*num_tri, 0)
    y_exponent = np.max(num_tri - 1, 0)
    return x**x_exponent * y**y_exponent
