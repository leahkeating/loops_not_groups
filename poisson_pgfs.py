# poisson pgfs

import numpy as np
import math

# We set the distribution p_st to two independent Poisson distribution
#   with means mu for single edges and nu for triangles
def G0(x, y, mu, nu):
    value = 0
    for s in range(25):
        for t in range(25):
            value+=(np.exp(-mu-nu)/(math.factorial(s)*math.factorial(t)))*(mu*x)**s*(nu*y)**t
    return value

def Gs(x, y, mu, nu):
    value = 0
    for s in range(25):
        for t in range(25):
            value+=s*(np.exp(-mu-nu)/(math.factorial(s)*math.factorial(t)))*(mu*x)**(s-1)*(nu*y)**t
    return value

def Gt(x, y, mu, nu):
    value = 0
    for s in range(25):
        for t in range(25):
            value+=t*(np.exp(-mu-nu)/(math.factorial(s)*math.factorial(t)))*(mu*x)**(s)*(nu*y)**(t-1)
    return value