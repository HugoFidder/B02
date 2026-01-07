import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from package4 import A_enclosed, y_linspace, t, b, zmax_front, zmax_rear, zmin_front, zmin_rear, Emod, x_arr, z_arr, A_arr, funcChord, I_xx, distributed_Ixx, M_func, V_func
from bendingtorsion import funcT
from constants import *

# shloks stuff
t_spar_root = .013
t_flage_root = .00039

#variables
k_s = 10
k_v = 2
k_c = 4
K = 3
pi = np.pi

""" Heights of front and rear spar """
h_f_spar = zmax_front - zmin_front
h_r_spar = zmax_rear - zmin_rear

# Divide by funcChord(0) to get relative thickness distribution, because t_spar_root are at root 
t_spar_arr = t_spar_root * funcChord(y_linspace)/funcChord(0)
t_flange_arr = t_flage_root * funcChord(y_linspace)/funcChord(0)

""" Stringer parameters """
a = .04    # height and width of L stringer section
t_s = .006   # thickness of L stringer section

""" Rib parameters """
#nribs = 5 #placeholder
#ribspacing = (b/2)/nribs
ribspacing = 0.092

""" Functions for L-stringer properties """
def L_area(a,t):
    return 2*a*t - t*t

def L_centroid(a,t):
    """ 
    Returns centroid measured from the corner of the L section. 
    The centroid is the same in x,y, because of symmetry.
    based on https://calcresource.com/moment-of-inertia-angle.html
    """
    return t * (a*t + a*a - t*t)/ (2*L_area(a,t))

def L_Imin(a,t):
    """ 
    Returns the smallest MoI for the L section, based on MoM book
    The taken axis the the axis which is normal to the symmetry axis
    """
    return t*a**3 / 12 + 2*(t*a)*(2*(.5*a+C_stringer)**2) 

#web buckling front spar

def b_plate_front(spanlocation): 
    for itvar in range(0,round(ribspacing/(b/2))+1):
        if spanlocation <= (y_linspace[-1]/(round(ribspacing/(b/2))+1))*(itvar+1):
            b_plate_f = min(ribspacing, h_f_spar*funcChord(((y_linspace[-1])/(round(ribspacing/(b/2))+1))*(itvar+1)))
    return b_plate_f

def b_plate_rear(spanlocation):
    for itvar in range(0,round(ribspacing/(b/2))+1):
        if spanlocation <= (y_linspace[-1]/(round(ribspacing/(b/2))+1))*(itvar+1):
            b_plate_r = min(ribspacing, h_r_spar*funcChord(((y_linspace[-1])/(round(ribspacing/(b/2))+1))*(itvar+1)))
    return b_plate_r

#front is a boolean designating whether its caluclating for front (true) or back (false) spar
def spar_buckling_stress(thickness, y_coord , front, k_s=10):
    """ 
    Critical shear stress for spar web buckling from manual 
    """
    if front == True:
        b_plate = b_plate_front(y_coord)
    else:
        b_plate = b_plate_rear(y_coord)
    return ((pi**2)*k_s*Emod)/(12*(1-poisson**2)) * ((thickness/b_plate)**2) #b was b_plate

def column_buckling_stress(K,I_min,A,L):
    """ 
    Euler buckling formula taken from manual
    """
    return (K*(pi**2)*Emod*I_min)/(L*L*A)

def skin_buckling_stress(t,y_coord, k_c=4):
    """ 
    Critical compressive stress for skin buckling from manual 
    """
    b = min(funcChord(y_coord), ribspacing) 
    return ((pi**2)*k_c*Emod)/(12*(1-poisson**2)) * ((t/b)**2)

# Interal forces
V_arr = V_func(y_linspace)
Torque_arr = np.empty_like(y_linspace)

# Torque integrated to get internal torque distribution
for i,y in enumerate(y_linspace):
    T = sp.integrate.quad(funcT,y,b/2)[0]
    Torque_arr[i] = T

c_arr = funcChord(y_linspace)

#second moment of inertia for I shape
C_stringer = L_centroid(a,t)
A_stringer = L_area(a,t_s)
I_min_stringer = L_Imin(a,t_s)

#spar web buckling
f_spar_crit_stress = []
for y in range(0, len(y_linspace)):
    f_spar_crit_stress.append(spar_buckling_stress(t_spar_arr[y], y_linspace[y], True))

r_spar_crit_stress = []
for y in range(0, len(y_linspace)):
    r_spar_crit_stress.append(spar_buckling_stress(t_spar_arr[y], y_linspace[y], False))

# f_spar_crit_stress = spar_buckling_stress(t_spar_arr, b_plate_front, k_s)
# # print(f'Critical stress for web of the front spar {tau_critical_front}')

# #web buckling rear spar
# r_spar_crit_stress = spar_buckling_stress(t_spar_arr, b_plate_rear, k_s)
# # print(f'Critical stress for web of the rear spar {tau_critical_rear}')

#average web shear  
tau_average = V_arr/((h_f_spar+h_r_spar) * funcChord(y_linspace)*t_spar_arr)
tau_max_average = k_v * tau_average

# shear due to torsion
q = Torque_arr/(2*A_enclosed)

tau_max_ave_front = tau_max_average + q/t_spar_arr
tau_max_ave_rear = tau_max_average - q/t_spar_arr

print("taumaxav", tau_max_average[3])
print("fsparstress", f_spar_crit_stress[3])
    
#skin buckling
sigma_crit_arr = []
for y in y_linspace:
    sigma_crit_arr.append(skin_buckling_stress(t, y))
# print(f'Critical stress for skin (top) {sigma_crit_top}')

nult = 3.11*1.5

sigma_max_compressive = nult*max(zmax_front,zmax_rear)*funcChord(y_linspace)*M_func(y_linspace) / distributed_Ixx(y_linspace)
# print(f'Applied stress compressive {sigma_max_compressive}')
sigma_max_tensile = min(zmin_front,zmin_rear)*funcChord(y_linspace)*M_func(y_linspace) / distributed_Ixx(y_linspace)
# print(f'Applied stress tensile {sigma_max_tensile}')

#column buckling
sigma_crit_buckling = column_buckling_stress(1/4,I_min_stringer,A_stringer,ribspacing)
print(f'Critical column buckling stress: {round(sigma_crit_buckling/1e6,2)} MPa')

#margins of safety (MoF) 
MoF_web_front = f_spar_crit_stress / tau_max_average
MoF_web_rear = r_spar_crit_stress / tau_max_average

MoF_skin_front_compressive = sigma_crit_arr/ sigma_max_compressive 
MoF_skin_front_tensile = sigma_crit_arr / sigma_max_tensile

MoF_stringer_compressive = sigma_crit_buckling / sigma_max_compressive
MoF_stringer_tensile = sigma_crit_buckling / sigma_max_tensile

print(f'Maximum stresses:\nFront spar: {round(tau_max_ave_front[0]/1e6,2)} MPa\nRear spar: {round(tau_max_ave_rear[0]/1e6,2)} MPa\nCompressive skin: {round(max(sigma_max_compressive)/1e6,2)} MPa\nTensile skin: {round(min(sigma_max_tensile)/1e6,2)} MPa')

print()

#diagrams
plt.close('all')
plt.plot(y_linspace, MoF_web_front, label='Margin of safety (front web)')
plt.plot(y_linspace, MoF_web_rear, label='Margin of safety (rear web)')
#plt.plot(y_linspace, MoF_skin_front_compressive, label='Margin of safety (skin)')
plt.ylim(-5000,5000)
plt.ylabel('Safety Margin')
plt.xlabel('Spanwise Location (m)')
plt.title('Safety Margin in Spar Web Distribution along Span')
plt.legend()
plt.show()
