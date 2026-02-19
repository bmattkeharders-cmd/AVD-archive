# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 12:59:18 2023

@author: bmatt
"""

import numpy as np
import matplotlib.pyplot as plt
import imageio

R = 287.05
gam = 1.4
cp = 1005

T0 = 220
P0 = 26*10**3
M0 = 0.85

g = 9.81

hf = 42800000
pi_c = 18
pi_cP = 1.8
Tt4 = 1450

n = 1000
beta = np.linspace(0.1,4,n)

def get_B(beta):
    B = beta / (beta + 1)
    return B

def get_tau_r(gam, M0):
    tau_r = 1+ ((gam-1)/2)*(M0**2)
    return tau_r

def get_tau_c(pi_c, gam):
    tau_c = pi_c**((gam-1)/gam)
    return tau_c

def get_tau_gam(Tt4, T0):
    tau_gam = Tt4 / T0
    return tau_gam

def get_tau_cP(gam, pi_cP):
    tau_cP = pi_cP**((gam-1)/gam)
    return tau_cP

def get_tau_t(beta, gam, M0, pi_cP, Tt4, T0):
    term1 = get_B(beta)*get_tau_cP(gam, pi_cP)
    term2 = (1-get_B(beta))*get_tau_c(pi_c, gam)
    num = term1 + term2 - 1
    dom = (1-get_B(beta))*get_tau_gam(Tt4, T0)
    tau_t = 1 - ((get_tau_r(gam, M0))*(num/dom))
    return tau_t

def get_M9overM0(beta, gam, M0, pi_cP, Tt4, T0):
    num = (get_tau_t(beta, gam, M0, pi_cP, Tt4, T0)*get_tau_r(gam, M0)*get_tau_c(pi_c, gam)) - 1
    dom = get_tau_r(gam, M0) - 1
    M9overM0 = np.sqrt(num / dom)
    return M9overM0

def get_T9overT0(Tt4, T0, gam, M0, pi_c):
    T9overT0 = (get_tau_gam(Tt4, T0)) / (get_tau_r(gam, M0)*get_tau_c(pi_c, gam))
    return T9overT0

def get_M9PoverM0(gam, M0, pi_cP):
    num = (get_tau_r(gam, M0) * get_tau_cP(gam, pi_cP)) -  1
    dom = get_tau_r(gam, M0) - 1
    M9PoverM0 = np.sqrt(num / dom)
    return M9PoverM0

def get_specThrust(gam, R, T0, M0, beta, pi_cP, Tt4, pi_c):
    LT = (np.sqrt(gam*R*T0)) * M0
    term1 = (1-get_B(beta))
    term2 = (get_M9overM0(beta, gam, M0, pi_cP, Tt4, T0)*get_T9overT0(Tt4, T0, gam, M0, pi_c)) - 1
    term3 = get_B(beta) * (get_M9PoverM0(gam, M0, pi_cP) - 1)
    specThrust = LT * ((term1*term2) + term3)
    return specThrust

specThrust = np.zeros(n)
for i in range(n):
    specThrust[i] = get_specThrust(gam, R, T0, M0, beta[i], pi_cP, Tt4, pi_c)

plt.plot(beta, specThrust)
plt.xlabel('Bypass ratio')
plt.ylabel('Specific Thrust, N/kg/s')
plt.xlim((0, np.max(beta)))
plt.ylim((0,1500))
plt.title("Specific Thrust vs Bypass Ratio")
plt.grid()

##############################################

def get_tau_f(hf, cp, T0):
    tau_f = (hf / (cp*T0))
    return tau_f

def get_f(beta, hf, cp, T0, Tt4, gam, M0, pi_c):
    num = 1- get_B(beta)
    dom = get_tau_f(hf, cp, T0) - get_tau_gam(Tt4, T0)
    TT = get_tau_gam(Tt4, T0) - (get_tau_r(gam, M0) * get_tau_c(pi_c, gam))
    f = (num/dom) * TT
    return f

def get_Isp(beta, hf, cp, T0, Tt4, gam, M0, pi_c, g, R, pi_cP):
    LE = get_specThrust(gam, R, T0, M0, beta, pi_cP, Tt4, pi_c)
    ME = 1 / g
    TE = 1 / get_f(beta, hf, cp, T0, Tt4, gam, M0, pi_c)
    Isp = LE*ME*TE
    return Isp

Isp = np.zeros(n)

for i in range(n):
    Isp[i] = get_Isp(beta[i], hf, cp, T0, Tt4, gam, M0, pi_c, g, R, pi_cP)


# def create_frame(t):
#     ax = plt.axes(projection = '3d')
#     ax.plot3D(beta, specThrust, Isp)
#     ax.view_init(20, t)
#     plt.savefig(f'.{t}.png', 
#                 transparent = False,  
#                 facecolor = 'white'
#                )
#     return

# angle = np.linspace(0,359,180)
# for t in angle:
#     create_frame(t)
    
# frames = []
# for t in angle:
#     image = imageio.imread(f'.{t}.png')
#     frames.append(image)

# imageio.mimsave('./thisOne.gif',
#                 frames,          
#                 fps = 30)

# plt.plot(beta, Isp)
# plt.xlabel('Bypass ratio')
# plt.ylabel('Specific Impulse, s')
# plt.xlim((0, np.max(beta)))
# plt.ylim((2000, 10000))
# plt.title("Specific Impulse vs Bypass Ratio")
# plt.grid()

