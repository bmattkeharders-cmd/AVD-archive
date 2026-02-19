# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 13:02:57 2023

@author: yohana
"""
import math
import numpy as np
import matplotlib.pyplot as plt

n = 1000
M = np.linspace(.01,.99,n)
beta = np.zeros(n)

def get_beta():
    for i in range(n):
        beta[i] = 1-M[i]**2

eta = np.zeros(n)
def get_eta():
    for i in range(n):
        eta[i] = 8.594/(2*3.14/beta[i])
    

AR = 1.94

gamma = 31
Sexp = 27
Sref = 27.5
d = 2
b = 10

get_beta()
get_eta()
CLa = np.zeros(n)
Cla_in = 0.5

def get_CLa():
    for i in range(n):
        CLa[i] = (2*3.14*AR)/((2+(4+((AR**2*beta[i]**2)/(eta[i]**2))*(1+(math.tan(gamma)**2)/(beta[i]**2)))**.5)) 
        # if .005 >= CLa[i] - CLa_in:
        #     idx = i
        # CLa_in = CLa[idx]
    
get_CLa()
plt.plot(M, CLa)
plt.grid()
plt.ylabel('Coefficient of Lift Slope')
plt.xlabel('Mach number')


