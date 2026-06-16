# Loops, not groups theory code for probability and size

import numpy as np
import negative_binomial_pgfs as pgfs
import os

# Run from this file's directory so the data/ subfolder resolves correctly
os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.makedirs("data", exist_ok=True)

###
# functions to load
###

def one_step(state, T, num_tri, deg, alpha):
    """Apply your map once and return (next_state, delta_rho)."""
    rho, u, v = state
    T2 = 1.0 - (1.0 - T)**alpha
    u2 = (1.0 - T) + T*Gs(u, v, num_tri, deg)
    v2 = (rho**2)*(1.0-T)*((1.0-T)*(1.0-T2)+T*(1.0-T))+2*rho*(1.0-rho)*(1.0-T)*((1.0-T)+T*(1.0-T))+(1.0-rho)**2
    rho2 = 1.0 - Gt(u2, v2, num_tri, deg)
    return np.array([rho2, u2, v2], dtype=float), (rho2 - rho)

def _v_of_rho(rho, T, num_tri, deg, alpha):
    """v-component of the map; it depends only on rho."""
    return one_step(np.array([rho, 0.0, 0.0]), T, num_tri, deg, alpha)[0][2]


def _u_of_v(v, T, num_tri, deg, alpha, maxit=200, tol=1e-13):
    """Inner fixed point u = (1 - T) + T*Gs(u, v) for a fixed v."""
    u = 1.0
    for _ in range(maxit):
        u_new = one_step(np.array([0.0, u, v]), T, num_tri, deg, alpha)[0][1]
        if abs(u_new - u) < tol:
            break
        u = u_new
    return u


def split_state_for_T(T, num_tri, deg, alpha, n_scan=400, maxit=80):
    """
    Giant-component size of the unstable (middle) branch at fixed T, or NaN if
    there is no bistability at this T. Finds the middle root of the scalar
    self-consistency H(rho) = (1 - Gt(u(rho), v(rho))) - rho, where v is slaved
    to rho and u to v.
    """
    def H(rho):
        v = _v_of_rho(rho, T, num_tri, deg, alpha)
        u = _u_of_v(v, T, num_tri, deg, alpha)
        return (1.0 - Gt(u, v, num_tri, deg)) - rho, u, v

    rhos = np.linspace(0.0, 1.0, n_scan)
    Hs = np.array([H(r)[0] for r in rhos])

    # the unstable fixed point is the negative -> positive crossing of H
    a = b = None
    for i in range(len(rhos) - 1):
        if Hs[i] < 0.0 <= Hs[i + 1]:
            a, b = rhos[i], rhos[i + 1]
            break
    if a is None:
        return np.nan

    for _ in range(maxit):
        m = 0.5 * (a + b)
        if H(m)[0] < 0.0:
            a = m
        else:
            b = m
    root = 0.5 * (a + b)
    _, u, v = H(root)
    return 1.0 - G0(u, v, num_tri, deg)

###

Gs = pgfs.Gs
Gt = pgfs.Gt
G0 = pgfs.G0

#Values of triangles per node to explore
num_tri_vec = [1.0, 2.0, 3.0]
deg = 6
alpha = 10.0

print("Calculating upper branch...")
for num_tri in num_tri_vec:
    print(f"num_tri={num_tri}, deg={deg}")
    Tlist = np.linspace(0.00,0.4,500)
    Size = []
    Prob = []
    u_in_vals = []
    v_in_vals = []
    # Initialise once; warm-start each T from the previous converged solution
    rho_in1 = 0.99
    rho_out1 = 0.99
    u_in1 = 0.01
    u_out1 = 0.01
    v_out1 = 0.01
    v_in1 = 0.01
    #Iterate over occupation probability
    for T in reversed(Tlist):
        T2 = 1-(1-T)**alpha
        T_loop = T2
        T_group = T
        u_out2 = (1-T) + T*Gs(u_out1,v_out1,num_tri,deg)
        u_in2 = (1-T) + T*Gs(u_in1,v_in1,num_tri,deg)
        v_in2 = rho_in1**2*(1-T)*((1-T)*(1-T_loop)+T*(1-T_group)) + (1-T)*(T*(1-T_group)+(1-T))*2*rho_in1*(1-rho_in1) + (1-rho_in1)**2
        v_out2 = 2*T*(1-T)*(1-T_group)*Gt(u_out1,v_out1,num_tri,deg)+(2*T*(1-T)*T_group+T**2)*(Gt(u_out1,v_out1,num_tri,deg)**2)+(1-T)**2
        rho_in2 = 1-Gt(u_in2,v_in2,num_tri,deg)
        rho_out2 = 1-Gt(u_out2,v_out2,num_tri,deg)
        while  abs(u_in2-u_in1) > 10**(-5) or abs(u_out2-u_out1) > 10**(-5) or abs(v_in2-v_in1) > 10**(-5) or abs(v_out2-v_out1) > 10**(-5) or abs(rho_in2-rho_in1) > 10**(-5) or abs(rho_out2-rho_out1) > 10**(-5):
            rho_in1 = rho_in2.copy()
            rho_out1 = rho_out2.copy()
            u_in1 = u_in2.copy()
            u_out1 = u_out2.copy()
            v_in1 = v_in2.copy()
            v_out1 = v_out2.copy()
            u_in2 = (1-T) + T*Gs(u_in1,v_in1,num_tri,deg)
            u_out2 = (1-T) + T*Gs(u_out1,v_out1,num_tri,deg)
            # need to edit the first term to include T_group
            v_in2 = rho_in1**2*(1-T)*((1-T)*(1-T_loop)+T*(1-T_group)) + (1-T)*(T*(1-T_group)+(1-T))*2*rho_in1*(1-rho_in1) + (1-rho_in1)**2
            v_out2 = 2*T*(1-T)*(1-T_group)*Gt(u_out1,v_out1,num_tri,deg)+(2*T*(1-T)*T_group+T**2)*(Gt(u_out1,v_out1,num_tri,deg)**2)+(1-T)**2
            rho_in2 = 1-Gt(u_in2,v_in2,num_tri,deg)
            rho_out2 = 1-Gt(u_out2,v_out2,num_tri,deg)
        Size.append(1-G0(u_in2,v_in2,num_tri,deg))
        Prob.append(1-G0(u_out2,v_out2,num_tri,deg))
        u_in_vals.append(u_in2)
        v_in_vals.append(v_in2)
        # carry converged solution forward as warm-start for next T
        u_in1 = u_in2
        u_out1 = u_out2
        v_in1 = v_in2
        v_out1 = v_out2
        rho_in1 = rho_in2
        rho_out1 = rho_out2
    Tlist = np.flip(Tlist)
    Size = [0 if abs(x) < 1e-3 else x for x in Size]
    Prob = (Prob)

    output_data = np.column_stack((Tlist, Size, u_in_vals, v_in_vals))
    supercrit = np.array([i for i in output_data if i[1]>=0.00001])
    np.savetxt("data/size_loops_no_groups_num_tri_{}_upper_branch.txt".format(int(num_tri)), supercrit[:, :2], fmt="%.6f")
    output_data = np.column_stack((Tlist, Prob))
    np.savetxt("data/prob_loops_no_groups_num_tri_{}.txt".format(int(num_tri)),
                    output_data, fmt="%.6f")
    
    print("Calculating lower branch...")
# for num_tri in num_tri_vec:
#     print(f"num_tri={num_tri}, deg={deg}")
#     Tlist = np.linspace(0.00,0.4,500)
    Size = []
    Prob = []
    u_in_vals = []
    v_in_vals = []
    # Initialise once; warm-start each T from the previous converged solution
    rho_in1 = 0.01
    rho_out1 = 0.01
    u_in1 = 0.99
    u_out1 = 0.99
    v_out1 = 0.99
    v_in1 = 0.99
    #Iterate over occupation probability
    for T in reversed(Tlist):
        T2 = 1-(1-T)**alpha
        T_loop = T2
        T_group = T
        u_out2 = (1-T) + T*Gs(u_out1,v_out1,num_tri,deg)
        u_in2 = (1-T) + T*Gs(u_in1,v_in1,num_tri,deg)
        v_in2 = rho_in1**2*(1-T)*((1-T)*(1-T_loop)+T*(1-T_group)) + (1-T)*(T*(1-T_group)+(1-T))*2*rho_in1*(1-rho_in1) + (1-rho_in1)**2
        v_out2 = 2*T*(1-T)*(1-T_group)*Gt(u_out1,v_out1,num_tri,deg)+(2*T*(1-T)*T_group+T**2)*(Gt(u_out1,v_out1,num_tri,deg)**2)+(1-T)**2
        rho_in2 = 1-Gt(u_in2,v_in2,num_tri,deg)
        rho_out2 = 1-Gt(u_out2,v_out2,num_tri,deg)
        while  abs(u_in2-u_in1) > 10**(-5) or abs(u_out2-u_out1) > 10**(-5) or abs(v_in2-v_in1) > 10**(-5) or abs(v_out2-v_out1) > 10**(-5) or abs(rho_in2-rho_in1) > 10**(-5) or abs(rho_out2-rho_out1) > 10**(-5):
            rho_in1 = rho_in2.copy()
            rho_out1 = rho_out2.copy()
            u_in1 = u_in2.copy()
            u_out1 = u_out2.copy()
            v_in1 = v_in2.copy()
            v_out1 = v_out2.copy()
            u_in2 = (1-T) + T*Gs(u_in1,v_in1,num_tri,deg)
            u_out2 = (1-T) + T*Gs(u_out1,v_out1,num_tri,deg)
            # need to edit the first term to include T_group
            v_in2 = rho_in1**2*(1-T)*((1-T)*(1-T_loop)+T*(1-T_group)) + (1-T)*(T*(1-T_group)+(1-T))*2*rho_in1*(1-rho_in1) + (1-rho_in1)**2
            v_out2 = 2*T*(1-T)*(1-T_group)*Gt(u_out1,v_out1,num_tri,deg)+(2*T*(1-T)*T_group+T**2)*(Gt(u_out1,v_out1,num_tri,deg)**2)+(1-T)**2
            rho_in2 = 1-Gt(u_in2,v_in2,num_tri,deg)
            rho_out2 = 1-Gt(u_out2,v_out2,num_tri,deg)
        Size.append(1-G0(u_in2,v_in2,num_tri,deg))
        Prob.append(1-G0(u_out2,v_out2,num_tri,deg))
        u_in_vals.append(u_in2)
        v_in_vals.append(v_in2)
        # carry converged solution forward as warm-start for next T
        u_in1 = u_in2
        u_out1 = u_out2
        v_in1 = v_in2
        v_out1 = v_out2
        rho_in1 = rho_in2
        rho_out1 = rho_out2
    Tlist = np.flip(Tlist)
    Size = [0 if abs(x) < 1e-3 else x for x in Size]
    Prob = (Prob)

    output_data = np.column_stack((Tlist, Size, u_in_vals, v_in_vals))
    subcrit = np.array([i for i in output_data if i[1]<0.00001])

    eps = 1e-5
    tol = 1e-6
    maxit = 200

    upper_branch_lhs = supercrit[-1]
    lower_branch_rhs = subcrit[-1]

    T_vals = np.linspace(upper_branch_lhs[0] + eps, lower_branch_rhs[0] - eps, 500)
    rho_boundary = []

    for T in T_vals:
        state_star = split_state_for_T(T, num_tri, deg, alpha)
        rho_boundary.append(state_star)  # unstable-branch size (NaN outside bistable region)
    
    np.savetxt("data/loops_no_groups_rho_boundary_num_tri_{}.txt".format(int(num_tri)), np.column_stack((T_vals, rho_boundary)), fmt="%.6f")

    # lower stable branch (Size 0) reaches the lower fold T_c, where it meets
    # the bottom of the unstable branch; emit it across [0, T_c] (the sweep
    # above only recorded it below the upper fold)
    rb = np.asarray(rho_boundary, dtype=float)
    T_c = np.nanmax(T_vals[np.isfinite(rb)]) if np.isfinite(rb).any() else supercrit[-1, 0]
    lower_T = np.sort(Tlist[Tlist <= T_c])
    np.savetxt("data/size_loops_no_groups_num_tri_{}_lower_branch.txt".format(int(num_tri)), np.column_stack((lower_T, np.zeros_like(lower_T))), fmt="%.6f")

    