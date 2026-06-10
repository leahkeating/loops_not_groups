# Loops, not groups theory code for probability and size

import numpy as np
import poisson_pgfs as pgfs

Gs = pgfs.Gs
Gt = pgfs.Gt
G0 = pgfs.G0

#Values of average connectivity to explore
mu_vec = [4.5, 2.5, 0.5]
nu_vec = [1.0, 2.0, 3.0]
alpha = 10.0

for mu,nu in zip(mu_vec, nu_vec):
    Tlist = np.linspace(0.00,0.4,500)
    Size = []
    Prob = []
    #Iterate over occupation probability
    for T in reversed(Tlist):
        print(f"T={T:.4f}, Size={1-G0(u_in2,v_in2,mu,nu):.4f}" if 'u_in2' in dir() else f"T={T:.4f}, starting")
        T2 = 1-(1-T)**alpha
        T_loop = T2
        T_group = T
        rho_in1 = 0.5
        rho_out1 = 0.5
        u_in1 = 0.5
        u_out1 = 0.5
        v_out1 = 0.5
        v_in1 = 0.5
        u_out2 = (1-T) + T*Gs(u_out1,v_out1,mu,nu)
        u_in2 = (1-T) + T*Gs(u_in1,v_in1,mu,nu)
        v_in2 = 2*(1-T)*(1-T_loop)*rho_in1**2 + (1-T)*(T*(1-T_group)+(1-T))*2*rho_in1*(1-rho_in1) + (1-rho_in1)**2
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
            v_in2 = 2*(1-T)*(1-T_loop)*rho_in1**2 + (1-T)*(T*(1-T_group)+(1-T))*2*rho_in1*(1-rho_in1) + (1-rho_in1)**2
            v_out2 = 2*T*(1-T)*(1-T_group)*Gt(u_out1,v_out1,mu,nu)+(2*T*(1-T)*T_group+T**2)*(Gt(u_out1,v_out1,mu,nu)**2)+(1-T)**2
            rho_in2 = 1-Gt(u_in2,v_in2,mu,nu)
            rho_out2 = 1-Gt(u_out2,v_out2,mu,nu)
        Size.append(1-G0(u_in2,v_in2,mu,nu))
        Prob.append(1-G0(u_out2,v_out2,mu,nu))
    Tlist = np.flip(Tlist)
    Size = [0 if abs(x) < 1e-3 else x for x in Size]
    jump = np.max(np.nonzero(Size))
    Size_up = list(Size[0:jump])
    Size_bot = list(Size[jump+1:])
    Prob = (Prob)


    output_data = np.column_stack((Tlist, Size))
    np.savetxt("data/size_loops_not_groups_nu_{}.txt".format(int(nu)),
                    output_data, fmt="%.6f")
    output_data = np.column_stack((Tlist, Prob))
    np.savetxt("data/prob_loops_not_groups_nu_{}.txt".format(int(nu)),
                    output_data, fmt="%.6f")