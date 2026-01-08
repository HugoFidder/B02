# Material properties
Emod = 72.4e9       #Pa
Gmod = 28e9         #Pa (shear modulus)
rho = 2780          #kg/m^3
sigma_y = 450e6     #Pa (yield strength)
y_ultimate = 485e6
poisson = .33      #-

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

