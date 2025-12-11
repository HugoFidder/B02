import numpy as np
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt
import pandas as pd
import math


# custom functions and constants
from constants import *
from functions import *
from package4 import *

ks=3
bPlate = 4     #<------- idk

shearCrit = (math.pi**2*k*Emod)/(12*(1-poission**2))*((t/bPlate)**2)
hf = 1
hr = 2
shearAvg = V/(hf*t+hr*t)

