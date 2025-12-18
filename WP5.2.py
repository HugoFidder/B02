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


stressArr = np.empty_like(y_linspace)

y_linspaceShort = np.linspace(0,(b/2)*(3/4),250)

maxStress = 450

yVals = [0,2.14,4.28,6.42,8.56,10.7]
stressVals = [128.28,111.97,87.16,53.91,17.96,0.42]

Stress = sp.interpolate.CubicSpline(yVals,stressVals)
for i,y in enumerate(y_linspaceShort):
    stressArr[i]=Stress(y)

margin = maxStress/stressArr

plt.plot(y_linspaceShort,margin)
plt.xlabel("Span Location[m]")
plt.ylabel("Margin of Safety[-]")
plt.show()
