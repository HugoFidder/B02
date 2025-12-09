import numpy as np

# Geometry

Cr = 3.49
taper = 0.36
b = 21.4
L = b / 2

E = 72.4e9
nu = 0.33
G = E / (2 * (1 + nu))

deflection_max = 0.075 * L


theta_limit_deg = 10.0
theta_limit = np.deg2rad(theta_limit_deg)


tau_allow = 295e6


# load polynomials

def chord(y):
    """Chord distribution"""
    return Cr * (1 + 2 * ((taper - 1) / b) * y)


def M_y(y):
    """Bending moment distribution"""
    return (7.08064596e-1 * y**5
            + 5.51003578     * y**4
            - 1.42780357e2   * y**3
            - 1.42417049e4   * y**2
            + 2.77884990e5   * y
            - 1.33986390e6)


def T_y(y):
    """Torque distribution"""
    return (3.48522513e-1 * y**5
            - 7.93608453    * y**4
            + 1.01585931e2  * y**3
            - 2.01119233e3  * y**2
            + 2.16967150e4  * y
            - 7.13525564e4)


# Geometry

zf_top = 0.0787
zf_bot = -0.04123125878
zr_top = 0.05161872904
zr_bot = -0.02161872904

# Spar height
h_front_factor = zf_top - zf_bot
h_rear_factor  = zr_top - zr_bot

# average vertical height
h_factor = 0.5 * (h_front_factor + h_rear_factor)


def box_dims(y):
    """
    Return b_box, h_box at a point in span
    """
    c = chord(y)
    b_box = 0.4 * c
    h_box = h_factor * c
    return b_box, h_box


def wingbox_Ixx(y, t):
    """
    Compute I_xx of the rectangular wingbox at point of span y, with wall thickness t.
    """
    c = chord(y)

    # actual z positions
    zf_t = zf_top * c
    zf_b = zf_bot * c
    zr_t = zr_top * c
    zr_b = zr_bot * c

    # Front spar
    h_front = zf_t - zf_b
    A_front = t * h_front
    zc_front = 0.5 * (zf_t + zf_b)

    Ixx_front_local = (t * h_front**3) / 12.0
    Ixx_front = Ixx_front_local + A_front * zc_front**2

    # rear spar
    h_rear = zr_t - zr_b
    A_rear = t * h_rear
    zc_rear = 0.5 * (zr_t + zr_b)

    Ixx_rear_local = (t * h_rear**3) / 12.0
    Ixx_rear = Ixx_rear_local + A_rear * zc_rear**2

    # top skin
    b_box, h_box = box_dims(y)

    z_top = 0.5 * (zf_t + zr_t)
    A_top = t * b_box
    # local ixx
    Ixx_top = A_top * z_top**2

    # bottom skin same thing
    z_bot = 0.5 * (zf_b + zr_b)
    A_bot = t * b_box
    Ixx_bot = A_bot * z_bot**2

    Ixx_total = Ixx_front + Ixx_rear + Ixx_top + Ixx_bot
    return Ixx_total


def wingbox_torsion_props(y, t):
    """
    For given point y and box wall thickness t (all 4 walls),
    return:
        A_m(y): midline area
        J(y): torsion constant
    """
    b_box, h_box = box_dims(y)
    A_m = b_box * h_box

    # lengths of each wall
    L_f = h_box
    L_r = h_box
    L_t = b_box
    L_b = b_box

    # all walls thickness = t
    denom = L_f / t + L_r / t + L_t / t + L_b / t
    J = 4 * A_m**2 / denom

    return A_m, J


# bending

def compute_deflection(N, dy=0.01):
    """
    integrate:
    w'' = M / (E Ixx(y))
    to get tip deflection at y = L.
    """
    w = 0.0
    theta = 0.0
    y = 0.0

    while y <= L:
        c = chord(y)
        t = N * c
        Ixx = wingbox_Ixx(y, t)
        M = M_y(y)

        wpp = M / (E * Ixx)

        theta += wpp * dy
        w     += theta * dy
        y     += dy

    return abs(w)


# torsion

def compute_torsion(N, n_steps=2001):
    """
    For given thickness law t(y) = N * chord(y),
    compute:
      - tip twist theta_tip
      - maximum shear stress tau_max in any wall
    assuming same t on all 4 walls.
    """
    ys = np.linspace(0.0, L, n_steps)
    thetap = []           # dtheta/dy
    tau_max_all = 0.0

    for y in ys:
        T = T_y(y)
        c = chord(y)
        t = N * c

        A_m, J = wingbox_torsion_props(y, t)

        # twist rate
        thetap.append(T / (G * J))

        # shear flow and stresses
        q = T / (2 * A_m)

        tau_here = abs(q / t)  # tau = q/t
        if tau_here > tau_max_all:
            tau_max_all = tau_here

    theta_tip = np.trapz(thetap, ys)
    return theta_tip, tau_max_all


# thickness factor across (for iteration)

def find_thickness_factor():

    N_low, N_high = 0.00000000000001, 0.5 # wrt chord

    for i in range(50):
        N_mid = 0.5 * (N_low + N_high)

        # Bending deflection
        d = compute_deflection(N_mid)

        # Torsion
        theta_tip, tau_max = compute_torsion(N_mid)

        # Check constraints
        violates = False
        if d > deflection_max:
            violates = True
        if abs(theta_tip) > theta_limit:
            violates = True
        if tau_max > tau_allow:
            violates = True

        if violates:
            # too thin -> increase thickness
            N_low = N_mid
        else:
            # safe -> try thinner
            N_high = N_mid

    # return conservative upper value
    return N_high


def thickness_distribution(N, n_pts=201):
    """
    Build span-wise list of thickness for t(y) = N * chord(y).
    """
    ys = np.linspace(0.0, L, n_pts)
    ts = [N * chord(y) for y in ys]
    return ys.tolist(), ts


# run

N_opt = find_thickness_factor()
print(f"Optimal thickness factor N = {N_opt:.6f} (t = N * chord)")


d_final = compute_deflection(N_opt)
theta_tip_final, tau_max_final = compute_torsion(N_opt)

print(f"Tip deflection = {d_final:.4f} m "
        f"({100*d_final/L:.2f}% of semi-span)")
print(f"Tip twist      = {np.rad2deg(theta_tip_final):.3f} deg")
print(f"Max shear tau  = {tau_max_final/1e6:.2f} MPa")

# Span-wise thickness distribution
y_list, t_list = thickness_distribution(N_opt, n_pts=201)


print(f"Root thickness t(0)   = {1e3 * t_list[0]:.2f} mm")
print(f"Tip thickness  t(L)   = {1e3 * t_list[-1]:.2f} mm")

