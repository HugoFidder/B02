import numpy as np
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt
import pandas as pd

# Material properties
Emod = 72.4e9       #Pa
Gmod = 28e9         #Pa (shear modulus)
rho = 2780          #kg/m^3
sigma_y = 450e6     #Pa (yield strength)
poission = .33      #-

# Environment
M = .68
MTOM = 10689
OEM = 6233.61
# T = 200
rhoCruise = 0.28736
V = 200.6
rhoISA = 1.225
qCruise = rhoCruise*(V**2)/2

# Wing properties
quarterSweep = 11       #deg
taper = 0.36            #-
b = 21.4                #m
Cr = 3.49               #m
Ct = 1.26               #m
dihedral = 4            #deg
MAC = 2.55              #m
Xlemac = 6.14           #m
Xtemac = 8.69           #m
MACspanwise = 4.51      #m
Sw = 50.88              #m^2
Swf = 24.05             #m^2
GearDist = 2.58         #m
GearLength = 1.3        #m
Mwing = 964.5           #kg
Mgear = 189.28/2        #kg
Fuel = 2                #m^3 To be determined
t = .005                #m/c (spar thickness)
x_front = 0.3           #-
x_rear = 0.7            #-
zmax_front = 0.0787        #found from NACA 2412 geometry at 0.3c
zmin_front = -0.04123125878      #-//-
zmax_rear = 0.05161872904        #found from NACA 2412 geometry at 0.7c
zmin_rear = -0.02161872904         #-//-
CL0 = 0.180034          #-
CL10 = 1.042209         #-
CM0 = -0.166919         #-
CM10 = -0.783209        #-
alpha_to = 21.05        #deg
alpha_land = 18.66      #deg
alpha_cruise = 0.91     #deg 
Mfuel = 3167            #kg
print ("WEIGHT SSDADSDDS",MTOM -Mfuel)

#Number of stringers on top and bottom
top_stringers = 10    
bottom_stringers = 10
N_stringers = top_stringers + bottom_stringers
I_xx = I_zz = I_xz = .0

#geometry inputs:
A_top = .00024
A_bottom = .0002

"""
x,z and area off all centroids of parts. One part has the same index number on the arrays.
I.e. the first stringer has x position x_arr[0], z position z_arr[0] and area A_arr[0]
The coordinates will be defined relative to the LE first and x,z will be relative to the chord (x=0 is LE and x=1 is TE)
"""

x_arr = np.zeros(N_stringers + 4,dtype=np.float64) # +4 accounts for the big spars
z_arr = np.zeros_like(x_arr) #zero_like generates an array of zeros with the same shape as x_arr
A_arr = np.zeros_like(x_arr)

def determine_CLalpha(CL10,CL0,dalpha):
    return (CL10 - CL0) / dalpha

def determine_CL(V,n,W=MTOM*9.81,rho=rhoISA,s=Sw):
    return n*W/(.5 * rho * V * V * s)

def determine_centroid(x_arr,z_arr,A_arr):
    # Calculate the centroid of the wingbox based on the geometric parameters and return as an array (x,z)
    X = sum(x_arr * A_arr)
    Z = sum(z_arr * A_arr)
    return (np.array((X,Z))/sum(A_arr))

def parallel_axis(area,coord1,coord2):
    # Calculate parallel-axis contribution based on areas and coords of spars and strigers
    return sum(area * coord1 * coord2)

def spar_length(coord1,coord2):
    return np.sqrt((coord2[0]-coord1[0])**2 + (coord2[1]-coord1[1])**2)

def midpoint(coord1,coord2):
    return [(coord1[0]+coord2[0])/2 , (coord1[1]+coord2[1])/2]

# Define coordinates of the 4 big spars
spars = np.empty(4, dtype=object) #Array containing the 4 coordinates of the saprs
spars[0] = np.array([x_front,zmin_front])  #Bottom left
spars[1] = np.array([x_front,zmax_front])   #Top left
spars[2] = np.array([x_rear,zmax_rear])    #Top right
spars[3] = np.array([x_rear,zmin_rear])   #Bottom right

 
# Calculate parameters for 4 big spars
for i in range(3):
    centroid = midpoint(spars[i],spars[i+1])
    L = spar_length(spars[i],spars[i+1])
    sin = (spars[i+1][1]-spars[i][1])/L
    cos = (spars[i+1][0]-spars[i][0])/L
    I_xx += t * L**3 / 12 * sin**2
    I_zz += t * L**3 / 12 * cos**2
    I_xz += t * L**3 / 12 * sin * cos
    A = t * L
    print(f"Point {i+1} to {i+2}: Centroid = {centroid}, Area = {A}, I_xx = {t * L**3 / 12 * sin**2}, I_zz = {t * L**3 / 12 * cos**2}, I_xz = {t * L**3 / 12 * sin * cos}")
    # Assign to arrays
    x_arr[N_stringers + i] = centroid[0]
    z_arr[N_stringers + i] = centroid[1]
    A_arr[N_stringers + i] = A

# Last spar from point 3 to point 0
centroid = midpoint(spars[3],spars[0])
L = spar_length(spars[3],spars[0])
sin = (spars[3][1]-spars[0][1])/L
cos = (spars[3][0]-spars[0][0])/L
I_xx += t * L**3 / 12 * sin**2
I_zz += t * L**3 / 12 * cos**2
I_xz += t * L**3 / 12 * sin * cos
A = t * L
print(f"Point {3} to {0}: Centroid = {centroid}, Area = {A}, I_xx = {t * L**3 / 12 * sin**2}, I_zz = {t * L**3 / 12 * cos**2}, I_xz = {t * L**3 / 12 * sin * cos}, sin = {sin}, cos = {cos}")

# Assign to arrays
x_arr[N_stringers + 3] = centroid[0]
z_arr[N_stringers + 3] = centroid[1]
A_arr[N_stringers + 3] = A

# Assign areas to area array corresponding to stringers on top and bottom spar
A_arr[0:top_stringers] = A_top                  # Top stringers                
A_arr[top_stringers:N_stringers] = A_bottom     # Bottom stringers

# Generate equidistant points to function as stringers
x_arr[0:top_stringers] = np.linspace(x_front,x_rear,top_stringers)
x_arr[top_stringers:N_stringers] = np.linspace(x_front,x_rear,top_stringers)

z_arr[0:top_stringers] = np.linspace(zmax_front,zmax_rear,top_stringers)
z_arr[top_stringers:N_stringers] = np.linspace(zmin_front,zmin_rear,bottom_stringers)

# Calulate centroid 
Cx,Cz = determine_centroid(x_arr,z_arr,A_arr)
print(f"The centroid is {Cx},{Cz}")
# Assume centroid is also Shear Centre for now
SCx,SCz = Cx,Cz

# Determine new coordinate system relative to centroid (xp = x' and zp = z')
xp_arr = x_arr - Cx 
zp_arr = z_arr - Cz

# Calculate final moments of inertia using parallel axis theorem
I_xx += parallel_axis(A_arr,zp_arr,zp_arr)
I_zz += parallel_axis(A_arr,xp_arr,xp_arr)
I_xz += parallel_axis(A_arr,xp_arr,zp_arr)

print(f"For the entire wingbox:\nI_xx = {I_xx}\nI_zz = {I_zz}\nI_xz={I_xz}")


# XFLR data from postive y onwards at aoa 0
df = pd.read_csv('./Alpha0.txt', on_bad_lines='skip')
Alpha0_arr = df.to_numpy()

y_span0 = Alpha0_arr[19:,0]
chord0 = Alpha0_arr[19:,1]
aInduced0 = Alpha0_arr[19:,2]
Cl0 = Alpha0_arr[19:,3]
Cdi0 = Alpha0_arr[19:,5]
Cm0 = Alpha0_arr[19:,7]

# XFLR data from postive y onwards at aoa 10
df1 = pd.read_csv('./Alpha10.txt', on_bad_lines='skip')
Alpha10_arr = df1.to_numpy()

y_span10 = Alpha10_arr[19:,0]
chord10 = Alpha10_arr[19:,1]
aInduced10 = Alpha10_arr[19:,2]
Cl10 = Alpha10_arr[19:,3]
Cdi10 = Alpha10_arr[19:,5]
Cm10 = Alpha10_arr[19:,7]

#Cl equation determination
Cl_intercept = CL0
Cl_slope = (Cl10-Cl0)/(CL10-CL0)

#Cm equation determination
Cm_slope = (Cm10-Cm0)/(CM10-CM0)
CM_a = (CM10-CM0)/10
CL_a = (CL10-CL0)/10


# Cl(y), CL, Cm(y), C(y) at alpha = 0
funcCl = sp.interpolate.CubicSpline(y_span0,Cl0)
funcCm = sp.interpolate.CubicSpline(y_span0,Cm0)

def funcChord(y):
    return Cr + (Ct - Cr)/(b/2) * y

funcCd = sp.interpolate.interp1d(y_span0, Cdi0,kind='cubic',fill_value="extrapolate")
# Extrapolated functions of dCl/dCL and dCm/dCM at every point.
funcdCldCL = sp.interpolate.CubicSpline(y_span0,Cl_slope)
funcdCmdCM = sp.interpolate.CubicSpline(y_span0,Cm_slope)
    

y_linspace = np.linspace(0,b/2,500)



#CM0+CM_a*(CL-CL0/CL_a)
def distributed_Cm(y,CM=CM0):
    return funcCm(y) + funcdCmdCM(y) * (CM - CM0)

def distributed_Cl(y,CL=CL0):
    return funcCl(y) + funcdCldCL(y) * (CL - CL0)


# Distributed second moments of area functions, scaled by local chord^4
# Because Ixx and Iyy are defined for chord = 1m

def distributed_Ixx(y):
    return I_xx * (funcChord(y))**4

def distributed_Izz(y):
    return I_zz * (funcChord(y))**4

def distributed_Ixz(y):
    return I_xz * (funcChord(y))**4


""" 
These function below should be changed to account for the changing local Cl as a function of big CL
(See above functions: distributed_Cm() and distributed_Cl())
"""

def lift(y,CL=CL0,q=qCruise):
    return distributed_Cl(y,CL) * funcChord(y) * q

Cd0 = 0.0488
def drag(y):
    return((funcCd(y)+Cd0/len(y_linspace))*qCruise*funcChord(y))



def moment(y,CL=CL0,CM=CM0, q=qCruise):
    print("ЕБАСИ",CM)
    c = funcChord(y)
    # Moment due to Cm
    Cm_m = distributed_Cm(y,CM) * (c**2) * q
    print(distributed_Cm(y,CM))
    # Moment due to Cl, assuming AC is at c/4
    Cl_m = distributed_Cl(y,CL) * c * q * (SCx - c/4)
    # print(Cm_m)
    return Cm_m + Cl_m


# Set up force arrays (V = shear, M = moment, D = dra)
V_arr = np.empty_like(y_linspace)
M_arr = np.empty_like(y_linspace)
D_arr = np.empty_like(y_linspace)
chordList = np.empty_like(y_linspace)
for i,y in enumerate(y_linspace):
    chordList[i] = funcChord(y)



#(integrated section volume/entire volue) * Mfuel
# volume equals area integrated over span
# area varies with the square of the chord
# chord varies linearly over the span
Mfuel/2

def Chord(y):
    return (Ct + (Cr-Ct)/(b/2)*(b/2-y))


def Area(y):
    return (Chord(y)**2) * sum(A_arr)

def IdilFuel(y):
   return((-funcChord(y)**2)*(-0.005007189*y+0.184802636)*775)

print("Fuel Weight",IdilFuel(b/2)-IdilFuel(0))
print(Mfuel/2)

#Critical Load Case
nFactor = -1 * 1.5 
vCrit = 231.4

def force_distribution(y_linspace,CL=CL0,q=qCruise,W_F_bool=False, nFactor =nFactor):
    L_arr = np.empty_like(y_linspace)
    Gear_W_array = np.empty_like(y_linspace)
    Fuel_W_array = np.empty_like(y_linspace)
    struct_W_array = np.empty_like(y_linspace)

    for i,y in enumerate(y_linspace):
        L_arr[i] = sp.integrate.quad(lift,y,b/2,args=(CL, q))[0]
        Fuel_W_array[i] = - (IdilFuel(b/2)-IdilFuel(y)) * 9.81*nFactor
        struct_W_array[i] = -sp.integrate.quad(Area, y, b/2)[0] * rho * 9.81*nFactor
        Gear_W_array[i] = - Mgear * 9.81*nFactor if (y<GearDist) else 0*nFactor
        #print("THE FUEL WEIGHT IS", F)
        # #D = sp.integrate.quad(drag,y,b/2)
        # A = sp.integrate.quad(Area, y, b/2)
        # # print(f"At {round(y / (b/2),2)}b is V = {L[0]}, M = {M[0]}, D = {D[0]}")
        # if (y<GearDist):
        #     V_arr[i] = L[0] + (Mgear+A[0]*rho+F) * 9.81
        # else:
        #     V_arr[i] = L[0] + (A[0]*rho+F)*9.81
        #     print("Lift", L)

    V_arr = L_arr + Gear_W_array + Fuel_W_array*W_F_bool+ struct_W_array

    # plt.plot(y_linspace,Fuel_W_array,label='Fuel Weight')
    # plt.plot(y_linspace,Gear_W_array,label='Gear Weight')
    # plt.plot(y_linspace,L_arr,label='Lift Force')
    plt.plot(y_linspace,V_arr,label='Total Shear Force')
    # plt.plot(y_linspace,struct_W_array,label='Structural Weight')
    plt.xlabel("Span (m)")
    plt.ylabel("Force (N)")
    plt.legend()    
    plt.title("Shear along Span")
    # plt.text(b/4, max(V_arr)*.9, f"CL = {round(CL,3)},q={round(q,3)}", fontsize=12)
    plt.show()
    V_arr_func = sp.interpolate.CubicSpline(y_linspace,V_arr)

    return V_arr_func,V_arr

    
    #D_arr[i] = D[0]




# #Plot the internal forces
# plt.subplot(211)
# plt.plot(y_linspace,V_arr)
# plt.subplot(212)
# plt.plot(y_linspace, M_arr)
# plt.show()
# plt.tight_layout()

#EOM to Find Internal Moment around x-axis

LiftTot = sp.integrate.quad(lift, 0, b/2) #integration of lift for one wing
def BendingMoment(x):
    return funcCl(x)*funcChord(x)*qCruise*x
TotalBendingMoment = sp.integrate.quad(BendingMoment,0,b/2)
print(TotalBendingMoment)


# funcV = sp.interpolate.CubicSpline(y_linspace,V_arr)

# BM_arr = np.empty_like(y_linspace)

# for i,y in enumerate(y_linspace):
#     BM_arr[0] = - TotalBendingMoment[0]
#     BM_arr[i] = sp.integrate.quad(funcV, y, b/2)[0] - TotalBendingMoment[0]

# print(BM_arr)

# Plot to verify interpolation functions
""" 
y_arr = np.linspace(0,b/2)
interp_Cm10 = np.zeros_like(y_arr)
interp_Cm5 = np.zeros_like(y_arr)
interp_Cm0 = np.zeros_like(y_arr)

interp_Cl10 = np.zeros_like(y_arr)
interp_Cl5 = np.zeros_like(y_arr)
interp_Cl0 = np.zeros_like(y_arr)

for i,y in enumerate(y_arr):
    interp_Cm10[i] = distributed_Cm(y,CM10)
    interp_Cm5[i] = distributed_Cm(y,(CM10 + CM0)/2)
    interp_Cm0[i] = distributed_Cm(y,CM0)
    interp_Cl10[i] = distributed_Cl(y,CL10)
    interp_Cl5[i] = distributed_Cl(y,(CL10 + CL0)/2)
    interp_Cl0[i] = distributed_Cl(y,CL0)

plt.subplot(211)
plt.plot(y_arr,interp_Cm10,label='interpolated cm10')
plt.plot(y_arr,interp_Cm5,label='interpolated cm5')
plt.plot(y_arr, interp_Cm0,label='interpolated cm0')
plt.plot(y_span0,Cm10,label='Cm10')
plt.plot(y_span0,Cm0,label='Cm0')
plt.title("Cm distribution verification")
plt.legend()


plt.subplot(212)
plt.plot(y_arr,interp_Cl10,label='interpolated cl10')
plt.plot(y_arr,interp_Cl5,label='interpolated cl5')
plt.plot(y_arr, interp_Cl0,label='interpolated cl0')
plt.plot(y_span0,Cl10,label='Cl10')
plt.plot(y_span0,Cl0,label='Cl0')
plt.title("Cl distribution verification")
plt.legend()

plt.show() 
"""

Load_CL = determine_CL(vCrit,nFactor, W = (MTOM-Mfuel)*9.81)
print(f"Critical Load Case CL: {Load_CL} at V = {vCrit} m/s and n = {nFactor}")
print("ALPHA",(Load_CL-0.180034)/((1.042209-0.180034)/10))

    


# Ixx_arr = distributed_Ixx(y_linspace)
# Izz_arr = distributed_Izz(y_linspace)
# Izx_arr = distributed_Ixz(y_linspace)
# plt.plot(y_linspace,Ixx_arr,label='Ixx')
# plt.plot(y_linspace,Izz_arr,label='Izz')
# plt.plot(y_linspace,Izx_arr,label='Izx')
# plt.xlabel("Span (m)")
# plt.ylabel("Moment of inertia (m^4)")
# plt.show()

V_func,V_arr = force_distribution(y_linspace,CL=Load_CL,q=.5*vCrit**2*rhoISA)

V_func2 = np.polyfit(y_linspace,V_arr,8)
print(f"The coefficients of the polynomial fit are: {V_func2}")

for i,y in enumerate(y_linspace):
    M0 = sp.integrate.quad(V_func,0,b/2)
    M = sp.integrate.quad(V_func,y,b/2)
    M_arr[i]=-M[0]

plt.plot(y_linspace,M_arr)
plt.title("Moment along Span")
plt.xlabel("Span (m)")
plt.ylabel("Moment (Nm)")
plt.show()

def deflection():
    # Correct integration: EI * d²v/dy² = M(y)
    # Integrate twice from root (y=0) with boundary conditions: v(0)=0, dv/dy(0)=0
    # First integration: dv/dy = ∫[0 to y] M(s)/(EI(s)) ds, with dv/dy(0) = 0
    curvature = M_arr / (Emod * distributed_Ixx(y_linspace))
    slope = sp.integrate.cumulative_trapezoid(curvature, y_linspace, initial=0)
    # Second integration: v = ∫[0 to y] slope(s) ds, with v(0) = 0
    v = sp.integrate.cumulative_trapezoid(slope, y_linspace, initial=0)
    return v

v = deflection()

print("Deflection at tip:", v[-1], "m")
plt.plot(y_linspace,v)
plt.title("Deflection along span")
plt.xlabel("Span (m)")
plt.ylabel("Deflection (m)")
# plt.gca().set_aspect('equal', adjustable='box')
plt.show()

A_enclosed = 0
for i in range(len(spars)-1):
    x0, z0 = spars[i]
    x1, z1 = spars[i+1]
    A_enclosed += (x0 * z1) - (x1 * z0)

A_enclosed = abs(A_enclosed) / 2
print("Enclosed area:", A_enclosed)
print(f"Enclosed area at cr = {Cr} m is {A_enclosed * Cr**2} m^2")

print("DEFLECTION: ",np.polyfit(y_linspace, v, 5))
print("SHEAR: ",np.polyfit(y_linspace, V_arr, 8))
print("MOMENT: ",np.polyfit(y_linspace, M_arr, 8))
