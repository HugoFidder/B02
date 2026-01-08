import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from package4 import chord0, lift, Cx, moment, Cz, c_arr, y_linspace, drag, b, Load_CL, CM0, CM_alpha, CL0, CL_alpha, vCrit, rhoISA, spars, Gmod, funcChord, t

# Ensure arrays
y_linspace = np.asarray(y_linspace)
c_arr = np.asarray(c_arr)

# Geometric positions
z_ac = 0.0171875 * c_arr
x_ac = 0.25 * c_arr
x_c  = Cx * c_arr
z_c  = Cz * c_arr

# Lift distribution
# L = np.asarray(lift(chordList))

# Drag distribution
D = np.asarray(drag(c_arr))


# Distance from centroid to aerodynamic centre
# hor_distance = x_ac - x_c
vert_dist = z_c - z_ac

# Moment from lift
# moment_ac_L = L * distance

# moment from drag
moment_ac_D = D*vert_dist
print(moment_ac_D)

#CM0+CM_a*(CL-CL0/CL_a)
# Cm pitching moment distribution
moment_ac_Cm = np.empty_like(y_linspace)
for i,y in enumerate(y_linspace):
    moment_ac_Cm[i] = np.asarray(moment(y, Load_CL, CM =CM0+CM_alpha*((Load_CL-CL0)/CL_alpha), q=.5*rhoISA*vCrit**2))

# Total torque distribution
moment_total = moment_ac_Cm + moment_ac_D

funcT = sp.interpolate.CubicSpline(y_linspace,moment_total)


Torque_arr = np.empty_like(y_linspace)
for i,y in enumerate(y_linspace):
    T = sp.integrate.quad(funcT,y,b/2)
    Torque_arr[i] = T[0]

# -------------------------
# Plot
# -------------------------
plt.figure()
plt.plot(y_linspace, Torque_arr, label="Total torque", linewidth=2)
plt.xlabel("Spanwise position y [m]")
plt.ylabel("Torque [Nm]")
plt.title("Spanwise Torque Distribution")
plt.legend()
plt.grid(True)
plt.tight_layout()
# plt.show()



Ttotal = sp.interpolate.CubicSpline(y_linspace,Torque_arr)
print()


def thetadz(y, Gmod = Gmod):
    return (1/(4*((0.038633743372*funcChord(y)**2)**2)*Gmod)*(0.994564938305/t)*funcChord(y))

def twist(y):
    return(thetadz(y)*Ttotal(y))

twist_arr = np.empty_like(y_linspace)
for i,y in enumerate(y_linspace):
    twist_arr[i] = sp.integrate.quad(twist, 0, y)[0]

plt.plot(y_linspace, twist_arr)
plt.ylabel("Twist [rad]")
plt.xlabel("Spanwise position y [m]")
plt.title("Wing twist over the Span")
# plt.show()


print("TORQUE: ",np.polyfit(y_linspace, Torque_arr, 5))
print("TWIST: ",np.polyfit(y_linspace, twist_arr, 5))


#dtheta/dz = (T/4Am**2*Gmod)*s/t          <--- to be determined after cross section is finalized

# Shoelace formula for area using spars list
