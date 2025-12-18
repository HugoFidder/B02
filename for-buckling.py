import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from package4 import chordList, A_enclosed, y_linspace, t, b, zmax_front, zmax_rear, zmin_front, zmin_rear, Emod, force_distribution, x_arr, z_arr, A_arr, funcChord, I_xx, distributed_Ixx, M_func
from bendingtorsion import funcT
from constants import *

#variables
k_s = 10
k_v = 2
k_c = 4
K = 3

# geometric spar parameters
a = .02
t = .004

def L_area(a,t):
    return 2*a*t - t*t

def L_centroid(a,t):
    C = t * (a*t + a*a - t*t)/ (2*L_area(a,t))
    return C

C = L_centroid(a,t)

def L_Imin(a,t):
    return t*a**3 / 12 + 2*(t*a)*(2*(.5*a+C)**2) 


#input:
nribs = 5 #placeholder
ribspacing = b/nribs

#forces
# V_arr_func, V_arr = force_distribution(y_linspace)
V_arr = V_func(y_linspace)

Torque_arr = np.empty_like(y_linspace)
for i,y in enumerate(y_linspace):
    T = sp.integrate.quad(funcT,y,b/2)
    Torque_arr[i] = T[0]

#second moment of inertia for I shape
A_spar = L_area(a,t)


#web buckling front spar
thickness = t * funcChord(y_linspace)
b_plate_front = min(ribspacing, (zmax_front-zmin_front)*funcChord(y_linspace))
tau_critical_front = ((3.14**2)*k_s*Emod)/(12*(1-poisson**2)) * ((thickness/b_plate_front)**2)

#web buckling rear spar
b_plate_rear = min(ribspacing, (zmax_rear-zmin_rear)*funcChord(y_linspace))
tau_critical_rear = ((3.14**2)*k_s*Emod)/(12*(1-poisson**2)) * ((thickness/b_plate_rear)**2)

#average web shear  
tau_average = V_arr/((zmax_front-zmin_front)*funcChord(y_linspace)*thickness + (zmax_rear-zmin_rear)*funcChord(y_linspace)*thickness)
tau_max_average = k_v * tau_average
q = Torque_arr/(2*A_enclosed)
tau_max_average += q

#skin buckling
sigma_crit_front = ((3.14**2)*k_c*Emod)/(12*(1-poisson**2))*((t/b_plate_front)**2)
sigma_crit_rear = ((3.14**2)*k_c*Emod)/(12*(1-poisson**2))*((t/b_plate_rear)**2)

sigma_max_skin = max(z_arr)*funcChord(y_linspace)*M_func(y_linspace) / distributed_Ixx(y_linspace)

#column buckling
I_stringer = 
sigma_crit_buckling = (K*(np.pi**2)*Emod*I)/((y_linspace**2)*A_spar)
sigma_max_buckling = 
# sigma_max 

# https://calcresource.com/moment-of-inertia-angle.html



I_min = L_Imin(a,t)

#column buckling
sigma_crit = (K*(3.12**2)*Emod*I_min)/((b/2)*A_spar)
print(f'Critical column buckling stress: {sigma_crit}')

#diagramms