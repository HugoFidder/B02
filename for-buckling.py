import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from package4 import chordList, A_enclosed, y_linspace, t, b, zmax_front, zmax_rear, zmin_front, zmin_rear, Emod, force_distribution
from bendingtorsion import funcT
# Ensure arrays
y_linspace = np.asarray(y_linspace)
chordList = np.asarray(chordList) #or funcChord?

#variables
k_s = 10
poisson = 0.3
b_plate_front = (zmax_front-zmin_front)
b_plate_rear = (zmax_rear-zmin_rear)
k_v = 2
k_c = 4
K = 3

#forces
V_arr_func, V_arr = force_distribution(y_linspace)

Torque_arr = np.empty_like(y_linspace)
for i,y in enumerate(y_linspace):
    T = sp.integrate.quad(funcT,y,b/2)
    Torque_arr[i] = T[0]

#second moment of inertia for I shape
I = 1*10**(-4)
A_spar = (10**(-2))


#web buckling front spar
tau_critical_front = ((3.14**2)*k_s*Emod)/(12*(1-poisson**2)) * ((t/b_plate_front)**2)

#web buckling rear spar
tau_critical_rear = ((3.14**2)*k_s*Emod)/(12*(1-poisson**2)) * ((t/b_plate_rear)**2)

#average web shear
tau_average = V_arr/(b_plate_front*t + b_plate_rear*t)
tau_max_average = k_v * tau_average
q = T/(2*A_enclosed)

#skin buckling
sigma_crit_front = ((3.14**2)*k_c*Emod)/(12*(1-poisson**2))*((t/b_plate_front)**2)
sigma_crit_rear = ((3.14**2)*k_c*Emod)/(12*(1-poisson**2))*((t/b_plate_rear)**2)

#column buckling
sigma_crit = (K*(3.12**2)*Emod*I)/((y_linspace**2)*A_spar)