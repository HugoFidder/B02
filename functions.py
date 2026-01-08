from constants import *
import numpy as np
import matplotlib.pyplot as plt

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
 
def spars_MoI(spars, t):

    """ 
    Calculate parameters for 4 big spars
    I_xx, I_zz, I_xz: Moments of inertia of the wingbox cross-section
    x_arr, z_arr, A_arr: Arrays containing x and z coordinates and area of each spar and stringer
    """

    I_xx = I_zz = I_xz = 0
    x_arr = np.zeros(4)
    z_arr = np.zeros(4)
    A_arr = np.zeros(4)
    
    for i0,i1 in ((0, 1), (1, 2), (2, 3), (3, 0)):
        centroid = midpoint(spars[i0],spars[i1])
        L = spar_length(spars[i0],spars[i1])
        sin = (spars[i1][1]-spars[i0][1])/L
        cos = (spars[i1][0]-spars[i0][0])/L
        I_xx += t * L**3 / 12 * sin**2
        I_zz += t * L**3 / 12 * cos**2
        I_xz += t * L**3 / 12 * sin * cos
        A = t * L
        # print(f"Point {i+1} to {i+2}: Centroid = {centroid}, Area = {A}, I_xx = {t * L**3 / 12 * sin**2}, I_zz = {t * L**3 / 12 * cos**2}, I_xz = {t * L**3 / 12 * sin * cos}")
        # Assign to arrays
        x_arr[i0] = centroid[0]
        z_arr[i0] = centroid[1]
        A_arr[i0] = A
        # print(f"x_arr: {x_arr}, z_arr: {z_arr}, A_arr: {A_arr}")
    return I_xx, I_zz, I_xz, x_arr, z_arr, A_arr
    