import numpy as np

taper = 0.36
Cr = 3.49
b = 21.4
L = b / 2
E = 72.4e9

dy = 0.01
deflection_max = 0.075 * L     # 7.5% of semi-span

# box geometry
ZF_TOP_FACTOR = 0.0787
ZF_BOT_FACTOR = -0.04123125878
ZR_TOP_FACTOR = 0.05161872904
ZR_BOT_FACTOR = -0.02161872904

# chordwise locations of spars as fraction of chord
X_FRONT_FACTOR = 0.3
X_REAR_FACTOR  = 0.7

# stringer area [m^2]
A_STR = 2*1.0e-4   # e.g. ~10 mm x 10 mm


# Moment
def Moment(y):
    return (7.08064596e-1*y**5 + 5.51003578*y**4
            - 1.42780357e2*y**3 - 1.42417049e4*y**2
            + 2.77884990e5*y - 1.33986390e6)


# chord
def chord(y):
    return Cr * (1 + 2 * ((taper - 1) / b) * y)


# ssection properties wingbox
def Calculator(c, t, top_stringers, bottom_stringers):
    """
    Compute section properties of a closed wingbox:

      - front spar (vertical)
      - rear spar (vertical)
      - top skin (horizontal)
      - bottom skin (horizontal)
      - optional top/bottom stringers

    """

    # z positions of spar ends
    zf_top = ZF_TOP_FACTOR * c
    zf_bot = ZF_BOT_FACTOR * c
    zr_top = ZR_TOP_FACTOR * c
    zr_bot = ZR_BOT_FACTOR * c

    # spar heights
    h_front = zf_top - zf_bot
    h_rear  = zr_top - zr_bot

    # chordwise positions of spars
    x_front = X_FRONT_FACTOR * c
    x_rear  = X_REAR_FACTOR * c


    b_box = (x_rear - x_front)          # 0.4 c
    z_top = 0.5 * (zf_top + zr_top)     # top skin z
    z_bot = 0.5 * (zf_bot + zr_bot)     # bottom skin z

    # spars

    # front spar
    A_front = t * h_front
    zc_front = 0.5 * (zf_top + zf_bot)
    xc_front = x_front

    Ixx_front_local = (t * h_front**3) / 12.0
    Izz_front_local = (h_front * t**3) / 12.0

    Ixx_front = Ixx_front_local + A_front * zc_front**2
    Izz_front = Izz_front_local + A_front * xc_front**2
    Ixz_front = A_front * xc_front * zc_front

    # rear spar
    A_rear = t * h_rear
    zc_rear = 0.5 * (zr_top + zr_bot)
    xc_rear = x_rear

    Ixx_rear_local = (t * h_rear**3) / 12.0
    Izz_rear_local = (h_rear * t**3) / 12.0

    Ixx_rear = Ixx_rear_local + A_rear * zc_rear**2
    Izz_rear = Izz_rear_local + A_rear * xc_rear**2
    Ixz_rear = A_rear * xc_rear * zc_rear


    # skins

    # top skin
    A_top = t * b_box

    Ixx_top = A_top * z_top**2

    xc_top = 0.5 * (x_front + x_rear)
    Izz_top = A_top * xc_top**2
    Ixz_top = A_top * xc_top * z_top

    # bottom skin
    A_bot = t * b_box
    Ixx_bot = A_bot * z_bot**2
    xc_bot = xc_top
    Izz_bot = A_bot * xc_bot**2
    Ixz_bot = A_bot * xc_bot * z_bot

    # sum

    I_xx = Ixx_front + Ixx_rear + Ixx_top + Ixx_bot
    I_zz = Izz_front + Izz_rear + Izz_top + Izz_bot
    I_xz = Ixz_front + Ixz_rear + Ixz_top + Ixz_bot

    # stringers

    # stringers along top skin
    if top_stringers > 0:
        # evenly spaced between spars, excluding at spar webs
        xs_top = np.linspace(x_front, x_rear, top_stringers + 2)[1:-1]
        for xs in xs_top:
            I_xx += A_STR * z_top**2
            I_zz += A_STR * xs**2
            I_xz += A_STR * xs * z_top

    # stringers along bottom skin
    if bottom_stringers > 0:
        xs_bot = np.linspace(x_front, x_rear, bottom_stringers + 2)[1:-1]
        for xs in xs_bot:
            I_xx += A_STR * z_bot**2
            I_zz += A_STR * xs**2
            I_xz += A_STR * xs * z_bot
    K = 0.0

    return I_xx, I_zz, I_xz, K


    # deflecgion solver 
def compute_deflection(N, top_stringers, bottom_stringers):

    w = 0.0          # deflection
    theta = 0.0      # slope

    y = 0.0
    while y <= L:

        c = chord(y)
        t = N * c

        I_xx, I_zz, I_xz, K = Calculator(c, t, top_stringers, bottom_stringers)

        M = Moment(y)
        wpp = M / (E * I_xx)

        theta += wpp * dy
        w     += theta * dy
        y += dy

    return abs(w)


# find thickness
def find_required_N(top_stringers):

    N_low  = 0.001   
    N_high = 0.15# 15% of chord

    for _ in range(40):

        N = 0.5 * (N_low + N_high)
        d = compute_deflection(N, top_stringers, top_stringers)

        if d > deflection_max:
            N_low = N # too flexible, increase thickness
        else:
            N_high = N # too stiff, reduce thickness

    return N



if __name__ == "__main__":
    thicknesses = []

    for n in range(31):

        N = find_required_N(n)

        t_root = N * Cr
        t_tip  = N * chord(L)

        thicknesses.append((n, t_root, t_tip))
        print(f"Stringers: {n}, "
              f"t_root = {1000*t_root:.2f} mm, "
              f"t_tip = {1000*t_tip:.2f} mm")

    print("\nSummary (stringers, t_root, t_tip):")
    for item in thicknesses:
        print(item)
