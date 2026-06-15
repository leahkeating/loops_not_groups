import os
import re

SRC = "/Users/leah/Documents/loops_not_groups"
DST = "/Users/leah/Documents/loops_not_groups/regular_networks"

# Robust replacement for the heuristic bisection used to trace the unstable
# (middle) branch. The original brackets a fixed pair of stale endpoints and
# bisects on the sign of a single map step; for the sharp regular-network
# transition that bracket collapses near the saddle node and returns NaN,
# leaving a gap so the dashed branch never reaches the upper stable solution.
#
# Here v is slaved to rho and u is slaved to v, so the 3D fixed point reduces to
# the scalar equation H(rho) = (1 - Gt(u(rho), v(rho))) - rho = 0. In the
# bistable region this has three roots; the unstable one is the middle root,
# where H crosses negative -> positive. We find it directly and return the
# giant-component size there (matching the upper/lower branch quantity), or NaN
# when there is no bistability at this T.
NEW_FUNCS = '''def _v_of_rho(rho, T, num_tri, deg, alpha):
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
'''

files = [
    "loops_no_groups_theory.py",
    "groups_no_loops_theory.py",
    "loops_and_groups_theory.py",
    "simple_contagion_theory.py",
]

import_block = (
    "import regular_pgfs as pgfs\n"
    "import os\n"
    "\n"
    "# Run from this file's directory so the data/ subfolder resolves correctly\n"
    "os.chdir(os.path.dirname(os.path.abspath(__file__)))\n"
    "os.makedirs(\"data\", exist_ok=True)"
)

for fname in files:
    with open(os.path.join(SRC, fname)) as f:
        text = f.read()

    # 1. swap the pgf import and make the script run from its own dir
    text = text.replace("import poisson_pgfs as pgfs", import_block)

    # 2. re-parametrize the sweep: two Poisson means -> num_tri + fixed deg
    text = text.replace(
        "#Values of average connectivity to explore",
        "#Values of triangles per node to explore",
    )
    text = text.replace(
        "mu_vec = [4.5, 2.5, 0.5]\nnu_vec = [1.0, 2.0, 3.0]",
        "num_tri_vec = [1.0, 2.0, 3.0]\ndeg = 6",
    )
    text = text.replace(
        "for mu,nu in zip(mu_vec, nu_vec):",
        "for num_tri in num_tri_vec:",
    )
    text = text.replace(
        'print(f"mu={mu}, nu={nu}")',
        'print(f"num_tri={num_tri}, deg={deg}")',
    )

    # 3. output filenames keyed by num_tri instead of nu
    text = text.replace("int(nu)", "int(num_tri)")
    text = text.replace("_nu_{}", "_num_tri_{}")

    # 4. all remaining mu/nu are PGF calls + helper signatures -> num_tri, deg
    text = text.replace("mu, nu", "num_tri, deg")
    text = text.replace("mu,nu", "num_tri,deg")

    # 5. drop dead plot-split lines (jump/Size_up/Size_bot are never read here,
    #    and np.max(np.nonzero(Size)) crashes on the all-subcritical lower branch
    #    that regular networks produce, since there is no bistable region)
    text = text.replace(
        "    jump = np.max(np.nonzero(Size))\n"
        "    Size_up = list(Size[0:jump])\n"
        "    Size_bot = list(Size[jump+1:])\n",
        "",
    )

    # 6. replace the heuristic unstable-branch bisection (stabilize_ic +
    #    split_state_for_T) with the robust scalar root-find above
    text = re.sub(
        r"def stabilize_ic\(.*?\n    return 0\.5\*\(hi \+ lo\)\n",
        NEW_FUNCS,
        text,
        flags=re.DOTALL,
    )

    # ...and update the boundary loop to the new signature, storing the unstable
    #    branch's giant-component size (NaN outside the bistable region)
    text = text.replace(
        "        state_star = split_state_for_T(T, upper_branch_lhs, lower_branch_rhs,\n"
        "                                    num_tri, deg, alpha, tol=tol, maxit=maxit, warmup=2)\n"
        "        rho_boundary.append(state_star[0])  # rho component (may be NaN if not bracketed)",
        "        state_star = split_state_for_T(T, num_tri, deg, alpha)\n"
        "        rho_boundary.append(state_star)  # unstable-branch size (NaN outside bistable region)",
    )

    # 7. The lower stable branch is the trivial Size-0 fixed point and exists up to
    #    the lower fold T_c, but the descending warm-started sweep only logs it
    #    below the upper fold T_sn (above T_sn the low-IC sweep rides the upper
    #    branch). Drop the early truncated save and re-emit the lower branch across
    #    [0, T_c] after the boundary, so it meets the bottom of the unstable branch.
    model = fname[: -len("_theory.py")]
    early_lower_save = (
        "    subcrit_plot = subcrit[subcrit[:, 0] <= supercrit[-1, 0]]\n"
        f'    np.savetxt("data/size_{model}_num_tri_{{}}_lower_branch.txt".format(int(num_tri)), subcrit_plot[:, :2], fmt="%.6f")\n'
    )
    text = text.replace(early_lower_save, "")

    boundary_save = (
        f'    np.savetxt("data/{model}_rho_boundary_num_tri_{{}}.txt".format(int(num_tri)), np.column_stack((T_vals, rho_boundary)), fmt="%.6f")\n'
    )
    boundary_save_plus_lower = boundary_save + (
        "\n"
        "    # lower stable branch (Size 0) reaches the lower fold T_c, where it meets\n"
        "    # the bottom of the unstable branch; emit it across [0, T_c] (the sweep\n"
        "    # above only recorded it below the upper fold)\n"
        "    rb = np.asarray(rho_boundary, dtype=float)\n"
        "    T_c = np.nanmax(T_vals[np.isfinite(rb)]) if np.isfinite(rb).any() else supercrit[-1, 0]\n"
        "    lower_T = np.sort(Tlist[Tlist <= T_c])\n"
        f'    np.savetxt("data/size_{model}_num_tri_{{}}_lower_branch.txt".format(int(num_tri)), np.column_stack((lower_T, np.zeros_like(lower_T))), fmt="%.6f")\n'
    )
    text = text.replace(boundary_save, boundary_save_plus_lower)

    out = os.path.join(DST, fname)
    with open(out, "w") as f:
        f.write(text)
    print("wrote", out)
