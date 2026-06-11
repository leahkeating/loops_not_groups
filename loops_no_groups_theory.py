# Loops, not groups theory code for probability and size

import numpy as np
import poisson_pgfs as pgfs

###
# functions to load
###

def one_step(state, T, mu, nu, alpha):
    """Apply your map once and return (next_state, delta_rho)."""
    rho, u, v = state
    T2 = 1.0 - (1.0 - T)**alpha
    u2 = (1.0 - T) + T*Gs(u, v, mu, nu)
    v2 = (rho**2)*(1.0-T)*((1.0-T)*(1.0-T2)+T*(1.0-T))+2*rho*(1.0-rho)*(1.0-T)*((1.0-T)+T*(1.0-T))+(1.0-rho)**2
    rho2 = 1.0 - Gt(u2, v2, mu, nu)
    return np.array([rho2, u2, v2], dtype=float), (rho2 - rho)

def stabilize_ic(state, T, mu, nu, alpha, n=2):
    """Optionally pre-iterate a couple times to 'stabilize' ICs like you did."""
    s = np.array(state, dtype=float)
    for _ in range(n):
        s, _ = one_step(s, T, mu, nu, alpha)
    return s

def split_state_for_T(T, upper_branch_lhs, lower_branch_rhs, mu, nu, alpha,
                      tol=1e-6, maxit=200, warmup=2):
    """
    For a fixed T, bisect between upper/lower initial states to find
    the state whose one-step delta_rho is ~ 0 (i.e., boundary between
    upward vs downward convergence).
    """
    # Upper/lower initial conditions (nudged like in your code)
    up0 = np.array([
        upper_branch_lhs[1] - eps,  # rho
        upper_branch_lhs[2] + eps,  # u
        upper_branch_lhs[3] + eps   # v
    ], dtype=float)
    lo0 = np.array([
        lower_branch_rhs[1] + eps,  # rho
        lower_branch_rhs[2] - eps,  # u
        lower_branch_rhs[3] - eps   # v
    ], dtype=float)

    # Optional warmup passes (your “stabilize the ICs” step)
    up = stabilize_ic(up0, T, mu, nu, alpha, n=warmup)
    lo = stabilize_ic(lo0, T, mu, nu, alpha, n=warmup)

    # Evaluate directions
    _, dup = one_step(up, T, mu, nu, alpha)
    _, dlo = one_step(lo, T, mu, nu, alpha)

    # Ensure opposite signs; if not, try swapping or return NaN to flag it
    if dup == 0:
        return up
    if dlo == 0:
        return lo
    if np.sign(dup) == np.sign(dlo):
        # No sign change — boundary may be outside bracket; return NaNs
        return np.array([np.nan, np.nan, np.nan])

    hi, dhi = (up, dup)
    lo, dlo = (lo, dlo)

    # True bisection on the STATE (component-wise midpoint), using sign of Δρ
    for _ in range(maxit):
        mid = 0.5*(hi + lo)
        _, dmid = one_step(mid, T, mu, nu, alpha)

        if dmid == 0 or np.max(np.abs(hi - lo)) < tol:
            return mid

        if np.sign(dmid) == np.sign(dhi):
            hi, dhi = mid, dmid
        else:
            lo, dlo = mid, dmid

    # If we exit by iteration cap, return best midpoint
    return 0.5*(hi + lo)

###

Gs = pgfs.Gs
Gt = pgfs.Gt
G0 = pgfs.G0

#Values of average connectivity to explore
mu_vec = [4.5, 2.5, 0.5]
nu_vec = [1.0, 2.0, 3.0]
alpha = 10.0

print("Calculating upper branch...")
for mu,nu in zip(mu_vec, nu_vec):
    print(f"mu={mu}, nu={nu}")
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
        u_out2 = (1-T) + T*Gs(u_out1,v_out1,mu,nu)
        u_in2 = (1-T) + T*Gs(u_in1,v_in1,mu,nu)
        v_in2 = rho_in1**2*(1-T)*((1-T)*(1-T_loop)+T*(1-T_group)) + (1-T)*(T*(1-T_group)+(1-T))*2*rho_in1*(1-rho_in1) + (1-rho_in1)**2
        v_out2 = 2*T*(1-T)*(1-T_group)*Gt(u_out1,v_out1,mu,nu)+(2*T*(1-T)*T_group+T**2)*(Gt(u_out1,v_out1,mu,nu)**2)+(1-T)**2
        rho_in2 = 1-Gt(u_in2,v_in2,mu,nu)
        rho_out2 = 1-Gt(u_out2,v_out2,mu,nu)
        while  abs(u_in2-u_in1) > 10**(-5) or abs(u_out2-u_out1) > 10**(-5) or abs(v_in2-v_in1) > 10**(-5) or abs(v_out2-v_out1) > 10**(-5) or abs(rho_in2-rho_in1) > 10**(-5) or abs(rho_out2-rho_out1) > 10**(-5):
            rho_in1 = rho_in2.copy()
            rho_out1 = rho_out2.copy()
            u_in1 = u_in2.copy()
            u_out1 = u_out2.copy()
            v_in1 = v_in2.copy()
            v_out1 = v_out2.copy()
            u_in2 = (1-T) + T*Gs(u_in1,v_in1,mu,nu)
            u_out2 = (1-T) + T*Gs(u_out1,v_out1,mu,nu)
            # need to edit the first term to include T_group
            v_in2 = rho_in1**2*(1-T)*((1-T)*(1-T_loop)+T*(1-T_group)) + (1-T)*(T*(1-T_group)+(1-T))*2*rho_in1*(1-rho_in1) + (1-rho_in1)**2
            v_out2 = 2*T*(1-T)*(1-T_group)*Gt(u_out1,v_out1,mu,nu)+(2*T*(1-T)*T_group+T**2)*(Gt(u_out1,v_out1,mu,nu)**2)+(1-T)**2
            rho_in2 = 1-Gt(u_in2,v_in2,mu,nu)
            rho_out2 = 1-Gt(u_out2,v_out2,mu,nu)
        Size.append(1-G0(u_in2,v_in2,mu,nu))
        Prob.append(1-G0(u_out2,v_out2,mu,nu))
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
    jump = np.max(np.nonzero(Size))
    Size_up = list(Size[0:jump])
    Size_bot = list(Size[jump+1:])
    Prob = (Prob)

    output_data = np.column_stack((Tlist, Size, u_in_vals, v_in_vals))
    supercrit = np.array([i for i in output_data if i[1]>=0.00001])
    np.savetxt("data/size_loops_no_groups_nu_{}_upper_branch.txt".format(int(nu)), supercrit[:, :2], fmt="%.6f")
    output_data = np.column_stack((Tlist, Prob))
    np.savetxt("data/prob_loops_no_groups_nu_{}.txt".format(int(nu)),
                    output_data, fmt="%.6f")
    
    print("Calculating lower branch...")
# for mu,nu in zip(mu_vec, nu_vec):
#     print(f"mu={mu}, nu={nu}")
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
        u_out2 = (1-T) + T*Gs(u_out1,v_out1,mu,nu)
        u_in2 = (1-T) + T*Gs(u_in1,v_in1,mu,nu)
        v_in2 = rho_in1**2*(1-T)*((1-T)*(1-T_loop)+T*(1-T_group)) + (1-T)*(T*(1-T_group)+(1-T))*2*rho_in1*(1-rho_in1) + (1-rho_in1)**2
        v_out2 = 2*T*(1-T)*(1-T_group)*Gt(u_out1,v_out1,mu,nu)+(2*T*(1-T)*T_group+T**2)*(Gt(u_out1,v_out1,mu,nu)**2)+(1-T)**2
        rho_in2 = 1-Gt(u_in2,v_in2,mu,nu)
        rho_out2 = 1-Gt(u_out2,v_out2,mu,nu)
        while  abs(u_in2-u_in1) > 10**(-5) or abs(u_out2-u_out1) > 10**(-5) or abs(v_in2-v_in1) > 10**(-5) or abs(v_out2-v_out1) > 10**(-5) or abs(rho_in2-rho_in1) > 10**(-5) or abs(rho_out2-rho_out1) > 10**(-5):
            rho_in1 = rho_in2.copy()
            rho_out1 = rho_out2.copy()
            u_in1 = u_in2.copy()
            u_out1 = u_out2.copy()
            v_in1 = v_in2.copy()
            v_out1 = v_out2.copy()
            u_in2 = (1-T) + T*Gs(u_in1,v_in1,mu,nu)
            u_out2 = (1-T) + T*Gs(u_out1,v_out1,mu,nu)
            # need to edit the first term to include T_group
            v_in2 = rho_in1**2*(1-T)*((1-T)*(1-T_loop)+T*(1-T_group)) + (1-T)*(T*(1-T_group)+(1-T))*2*rho_in1*(1-rho_in1) + (1-rho_in1)**2
            v_out2 = 2*T*(1-T)*(1-T_group)*Gt(u_out1,v_out1,mu,nu)+(2*T*(1-T)*T_group+T**2)*(Gt(u_out1,v_out1,mu,nu)**2)+(1-T)**2
            rho_in2 = 1-Gt(u_in2,v_in2,mu,nu)
            rho_out2 = 1-Gt(u_out2,v_out2,mu,nu)
        Size.append(1-G0(u_in2,v_in2,mu,nu))
        Prob.append(1-G0(u_out2,v_out2,mu,nu))
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
    jump = np.max(np.nonzero(Size))
    Size_up = list(Size[0:jump])
    Size_bot = list(Size[jump+1:])
    Prob = (Prob)

    output_data = np.column_stack((Tlist, Size, u_in_vals, v_in_vals))
    subcrit = np.array([i for i in output_data if i[1]<0.00001])
    subcrit_plot = subcrit[subcrit[:, 0] <= supercrit[-1, 0]]
    np.savetxt("data/size_loops_no_groups_nu_{}_lower_branch.txt".format(int(nu)), subcrit_plot[:, :2], fmt="%.6f")

    eps = 1e-5
    tol = 1e-6
    maxit = 200

    upper_branch_lhs = supercrit[-1]
    lower_branch_rhs = subcrit[-1]

    T_vals = np.linspace(upper_branch_lhs[0] + eps, lower_branch_rhs[0] - eps, 500)
    rho_boundary = []

    for T in T_vals:
        state_star = split_state_for_T(T, upper_branch_lhs, lower_branch_rhs,
                                    mu, nu, alpha, tol=tol, maxit=maxit, warmup=2)
        rho_boundary.append(state_star[0])  # rho component (may be NaN if not bracketed)
    
    np.savetxt("data/loops_no_groups_rho_boundary_nu_{}.txt".format(int(nu)), np.column_stack((T_vals, rho_boundary)), fmt="%.6f")

    