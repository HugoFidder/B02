import numpy as np
from scipy.interpolate import interp1d

Cr = 3.49
taper = 0.36
b = 21.4
L = b / 2

E = 72.4e9
nu = 0.33
G = E / (2 * (1 + nu))

# Bending deflection limit
deflection_max = 0.075 * L   # 7.5% of semi-span

# Torsion limits
theta_limit_deg = 10.0
theta_limit = np.deg2rad(theta_limit_deg)

tau_allow = 295e6 / 1.5  # with safety factor

# Buckling settings
SF_BUCK = 1.5
K_SKIN = 7.25  # plate buckling coefficient

# Integration step
DY = 0.01

def chord(y):
    """Chord distribution."""
    return Cr * (1 + 2 * ((taper - 1) / b) * y)


def M_y(y):
    """Bending moment distribution."""
    return (7.08064596e-1 * y**5
            + 5.51003578     * y**4
            - 1.42780357e2   * y**3
            - 1.42417049e4   * y**2
            + 2.77884990e5   * y
            - 1.33986390e6)


def T_y(y):
    """Torque distribution."""
    return (3.48522513e-1 * y**5
            - 7.93608453    * y**4
            + 1.01585931e2  * y**3
            - 2.01119233e3  * y**2
            + 2.16967150e4  * y
            - 7.13525564e4)


# wingbox geometry

# z-location factors
zf_top = 0.0787
zf_bot = -0.04123125878
zr_top = 0.05161872904
zr_bot = -0.02161872904

# spar height factors
h_front_factor = zf_top - zf_bot
h_rear_factor  = zr_top - zr_bot

# average vertical height factor
h_factor = 0.5 * (h_front_factor + h_rear_factor)


def box_dims(y):
    """
    Return (b_box, h_box) at span location y.

    b_box: horizontal distance between spars (top/bottom skin length)
    h_box: average vertical distance between skins
    """
    c = chord(y)
    b_box = 0.4 * c
    h_box = h_factor * c
    return b_box, h_box


# stringers

X_FRONT_FACTOR = 0.3
X_REAR_FACTOR  = 0.7

t_st   = 6.0e-3
a_leg  = 40.0e-3
b_leg  = 40.0e-3

A_STR = t_st * a_leg + t_st * b_leg - t_st * t_st
Ixx_STR_LOCAL = 3.142e-10


# 5 day bender (bending)

def Ixx_flange_and_stringers(y, t_flange, top_stringers, bottom_stringers):
    c = chord(y)

    zf_t = zf_top * c
    zf_b = zf_bot * c
    zr_t = zr_top * c
    zr_b = zr_bot * c

    b_box, _ = box_dims(y)

    x_front = X_FRONT_FACTOR * c
    x_rear  = X_REAR_FACTOR  * c

    z_top = 0.5 * (zf_t + zr_t)
    z_bot = 0.5 * (zf_b + zr_b)

    A_top = t_flange * b_box
    A_bot = t_flange * b_box

    Ixx_total = A_top * z_top**2 + A_bot * z_bot**2

    if top_stringers > 0:
        xs_top = np.linspace(x_front, x_rear, top_stringers + 2)[1:-1]
        for _ in xs_top:
            Ixx_total += Ixx_STR_LOCAL + A_STR * z_top**2

    if bottom_stringers > 0:
        xs_bot = np.linspace(x_front, x_rear, bottom_stringers + 2)[1:-1]
        for _ in xs_bot:
            Ixx_total += Ixx_STR_LOCAL + A_STR * z_bot**2

    return Ixx_total


def compute_deflection(N_flange, n_top, n_bot, dy=DY):
    w = 0.0
    theta = 0.0
    y = 0.0

    while y <= L:
        c = chord(y)
        t_flange = N_flange * c
        Ixx = Ixx_flange_and_stringers(y, t_flange, n_top, n_bot)
        M = M_y(y)

        wpp = M / (E * Ixx)
        theta += wpp * dy
        w     += theta * dy
        y     += dy

    return abs(w)

# one, two, buckle my shoe (buckling)

def sigma_cr_skin_plate(t, b_panel, k=K_SKIN):

    if b_panel <= 0 or t <= 0:
        return 0.0
    return (k * np.pi**2 * E) / (12.0 * (1.0 - nu**2)) * (t / b_panel)**2


def max_top_skin_buckling_utilisation(N_flange, n_top, n_points=2001):

    ys = np.linspace(0.0, L, n_points)
    util_max = 0.0

    for y in ys:
        c = chord(y)
        t_fl = N_flange * c

        b_box, _ = box_dims(y)
        b_panel = b_box / (n_top + 1)  # panel width between stringers

        Ixx = Ixx_flange_and_stringers(y, t_fl, n_top, n_top) 
        M = M_y(y)

        zf_t = zf_top * c
        zr_t = zr_top * c
        z_top = 0.5 * (zf_t + zr_t)

        sigma_top = M * z_top / Ixx 
        sigma_comp = abs(min(sigma_top, 0.0))

        sigma_cr = sigma_cr_skin_plate(t_fl, b_panel, k=K_SKIN) / SF_BUCK

        util = sigma_comp / sigma_cr if sigma_cr > 0 else 1e9
        util_max = max(util_max, util)

    return util_max

# torsion

def wingbox_torsion_props(y, t_flange, t_web):
    b_box, h_box = box_dims(y)
    A_m = b_box * h_box

    L_f = h_box
    L_r = h_box
    L_t = b_box
    L_b = b_box

    denom = (L_f / t_web +
             L_r / t_web +
             L_t / t_flange +
             L_b / t_flange)

    J = 4.0 * A_m**2 / denom
    return A_m, J


def compute_torsion(N_flange, N_web, n_steps=2001):
    ys = np.linspace(0.0, L, n_steps)
    thetap = []
    tau_web_max = 0.0

    for y in ys:
        T = T_y(y)
        c = chord(y)
        t_flange = N_flange * c
        t_web    = N_web    * c

        A_m, J = wingbox_torsion_props(y, t_flange, t_web)

        thetap.append(T / (G * J))

        q = T / (2.0 * A_m)
        tau_web_here = abs(q / t_web)
        tau_web_max = max(tau_web_max, tau_web_here)

    theta_tip = np.trapz(thetap, ys)
    return theta_tip, tau_web_max


# looooooop de loooop

def find_flange_factor(n_top, n_bot):
    """
    Find N_flange such that:
      - deflection <= deflection_max
      - top skin buckling utilisation <= 1
    """
    N_low, N_high = 1e-6, 0.5

    for _ in range(70):
        N_mid = 0.5 * (N_low + N_high)

        d = compute_deflection(N_mid, n_top, n_bot)
        util_buck = max_top_skin_buckling_utilisation(N_mid, n_top)

        violates = (d > deflection_max) or (util_buck > 1.0)

        if violates:
            N_low = N_mid
        else:
            N_high = N_mid

    return N_high


def find_web_factor_given_flange(N_flange):
    N_low, N_high = 1e-6, 0.5

    for _ in range(60):
        N_mid = 0.5 * (N_low + N_high)
        theta_tip, tau_web_max = compute_torsion(N_flange, N_mid)

        violates = (abs(theta_tip) > theta_limit) or (tau_web_max > tau_allow)

        if violates:
            N_low = N_mid
        else:
            N_high = N_mid

    return N_high


def thickness_distributions(N_flange, N_web, n_pts=201):
    ys = np.linspace(0.0, L, n_pts)
    t_flange = [N_flange * chord(y) for y in ys]
    t_web    = [N_web    * chord(y) for y in ys]
    return ys.tolist(), t_flange, t_web


# stress (literally me)

def build_normal_stress_functions(N_flange, n_top, n_bot, n_points=2001):
    ys = np.linspace(0.0, L, n_points)

    sigma_top_vals = []
    sigma_bot_vals = []

    for y in ys:
        c = chord(y)
        t_fl = N_flange * c

        Ixx = Ixx_flange_and_stringers(y, t_fl, n_top, n_bot)
        M = M_y(y)

        zf_t = zf_top * c
        zr_t = zr_top * c
        zf_b = zf_bot * c
        zr_b = zr_bot * c

        z_top = 0.5 * (zf_t + zr_t)
        z_bot = 0.5 * (zf_b + zr_b)

        sigma_top_vals.append(M * z_top / Ixx)
        sigma_bot_vals.append(M * z_bot / Ixx)

    sigma_top_fun = interp1d(ys, sigma_top_vals, kind="cubic", fill_value="extrapolate")
    sigma_bottom_fun = interp1d(ys, sigma_bot_vals, kind="cubic", fill_value="extrapolate")
    return sigma_top_fun, sigma_bottom_fun


# mass

def estimate_wingbox_mass(N_flange, N_web, n_top, n_bot, rho=2780.0, n_steps=2001):
    ys = np.linspace(0.0, L, n_steps)

    mass = 0.0
    for i in range(len(ys) - 1):
        ym = 0.5 * (ys[i] + ys[i + 1])
        dy = ys[i + 1] - ys[i]

        c = chord(ym)
        b_box, h_box = box_dims(ym)

        t_fl = N_flange * c
        t_w  = N_web    * c

        A_skins = 2.0 * b_box * t_fl
        A_webs  = 2.0 * h_box * t_w
        A_strs  = (n_top + n_bot) * A_STR

        mass += rho * (A_skins + A_webs + A_strs) * dy

    return mass


# inter rivet (why did we choose to do this)

def rivet_pitch_from_inter_rivet_buckling(
    N_flange,
    n_top,
    n_bot,
    k_ir=1.0,
    SF_ir=1.5,
    n_points=2001,
    sigma_floor=1e3
):

    ys = np.linspace(0.0, L, n_points)
    s_max = np.zeros_like(ys)
    sigma_c = np.zeros_like(ys)

    C = (k_ir * np.pi**2 * E) / (12.0 * (1.0 - nu**2))

    for i, y in enumerate(ys):
        c = chord(y)
        t = N_flange * c

        Ixx = Ixx_flange_and_stringers(y, t, n_top, n_bot)
        M = M_y(y)

        z_top = 0.5 * (zf_top * c + zr_top * c)

        sigma_top = M * z_top / Ixx
        sigma_comp = abs(min(sigma_top, 0.0))
        sigma_use = max(sigma_comp, sigma_floor)

        sigma_c[i] = sigma_comp
        s_max[i] = t * np.sqrt(C / (SF_ir * sigma_use))

    s_design = float(np.min(s_max))
    y_crit = float(ys[np.argmin(s_max)])

    return s_design, y_crit, ys, s_max, sigma_c



# run loop

if __name__ == "__main__":

    TOP_RANGE = range(0, 21)
    BOT_RANGE = range(0, 21)
    ENFORCE_SYMMETRIC = True

    best = None

    print("Running sizing sweep over stringer counts (with top-skin buckling)...\n")
    print(f"Deflection limit: {deflection_max:.3f} m ({100*deflection_max/L:.2f}% of semi-span)")
    print(f"Twist limit     : {theta_limit_deg:.2f} deg")
    print(f"Web tau allow   : {tau_allow/1e6:.2f} MPa")
    print(f"Skin buckling   : k={K_SKIN:.2f}, SF={SF_BUCK:.2f}\n")

    for n_top in TOP_RANGE:
        bot_iter = [n_top] if ENFORCE_SYMMETRIC else BOT_RANGE
        for n_bot in bot_iter:

            N_flange = find_flange_factor(n_top, n_bot)
            N_web    = find_web_factor_given_flange(N_flange)

            d_final = compute_deflection(N_flange, n_top, n_bot)
            theta_tip_final, tau_web_final = compute_torsion(N_flange, N_web)
            util_buck = max_top_skin_buckling_utilisation(N_flange, n_top)

            ok = (
                (d_final <= deflection_max)
                and (abs(theta_tip_final) <= theta_limit)
                and (tau_web_final <= tau_allow)
                and (util_buck <= 1.0)
            )
            if not ok:
                continue

            mass = estimate_wingbox_mass(N_flange, N_web, n_top, n_bot)

            cand = {
                "n_top": n_top,
                "n_bot": n_bot,
                "N_flange": N_flange,
                "N_web": N_web,
                "deflection": d_final,
                "twist": theta_tip_final,
                "tau_web": tau_web_final,
                "buck_util": util_buck,
                "mass": mass
            }

            if (best is None) or (cand["mass"] < best["mass"]):
                best = cand

    if best is None:
        print("No feasible design found in the search range.")
        raise SystemExit

    n_top = best["n_top"]
    n_bot = best["n_bot"]
    N_flange = best["N_flange"]
    N_web = best["N_web"]

    print("--- BEST DESIGN FOUND ---")
    print(f"Top stringers    = {n_top}")
    print(f"Bottom stringers = {n_bot}")
    print("Stringer type: 40 x 40 x 6 mm L\n")

    print("\n--- Final performance ---")
    print(f"Tip deflection = {best['deflection']:.4f} m ({100*best['deflection']/L:.2f}% of semi-span)")
    print(f"Tip twist      = {np.rad2deg(best['twist']):.3f} deg")
    print(f"Max web shear τ= {best['tau_web']/1e6:.2f} MPa")
    print(f"Estimated semi-span wingbox mass = {best['mass']:.2f} kg\n")

    y_list, t_flange_list, t_web_list = thickness_distributions(N_flange, N_web, n_pts=201)

    print("--- Root & Tip thicknesses ---")
    print(f"Flange root t_f(0) = {1e3 * t_flange_list[0]:.2f} mm")
    print(f"Flange tip  t_f(L) = {1e3 * t_flange_list[-1]:.2f} mm")
    print(f"Web    root t_w(0) = {1e3 * t_web_list[0]:.2f} mm")
    print(f"Web    tip  t_w(L) = {1e3 * t_web_list[-1]:.2f} mm")

    sigma_top_fun, sigma_bottom_fun = build_normal_stress_functions(N_flange, n_top, n_bot)

    sample_ys = np.linspace(0.0, L, 6)
    print("\n--- Normal stress distribution (bending only) ---")
    for y_val in sample_ys:
        sig_top = float(sigma_top_fun(y_val)) / 1e6
        sig_bot = float(sigma_bottom_fun(y_val)) / 1e6
        print(f"y = {y_val:6.2f} m : σ_top = {sig_top:8.2f} MPa (top), σ_bottom = {sig_bot:8.2f} MPa (bottom)")

s_design, y_crit, ys, smax, sigc = rivet_pitch_from_inter_rivet_buckling(
    N_flange,
    n_top=best["n_top"],
    n_bot=best["n_bot"],
    k_ir=1.0,
    SF_ir=1.5
)

print("")
print(f"Inter-rivet buckling pitch limit: s <= {1000*s_design:.1f} mm")
