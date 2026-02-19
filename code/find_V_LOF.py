# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 19:31:45 2023

@author: bmatt
"""

import numpy as np
import matplotlib.pyplot as plt

cl_alpha = 2*np.pi
a = 340
AR1 = 5
AR2 = 3.5
n = 1000
V_arr = np.linspace(40,100,n)

cL0 = 0
tilt_angle = 8
weight = 6500*9.81
rho = 1.225
Sw1 = 20 

def get_cl(AR, cL0, V, angle):
    term1 = cl_alpha / (cl_alpha/(np.pi*AR))
    term2 = 1 / np.sqrt((1-(V**2/a**2)))
    term3 = angle*np.pi / 180
    cl = term1*term2*term3 + cL0
    return cl

def get_cl_needed(Sw, W, V):
    cl_need = (2*W) / (rho*Sw*V**2)
    return cl_need

cl = np.zeros(n)
for i in range(n):
    cl[i] = get_cl(AR1, cL0, V_arr[i], tilt_angle)

cl_need = np.zeros(n)
for i in range(n):
    cl_need[i] = get_cl_needed(Sw1, weight, V_arr[i])
    
plt.plot(V_arr, cl)
plt.plot(V_arr, cl_need)
plt.grid()

### intercept of these graphs is coefficient of list and speed at which lift off is possible. ####
## avoid tilt angles over 10 cause that's where stall occurs ###

## I used this to calculate approx take off distances. Pretty accurate #####