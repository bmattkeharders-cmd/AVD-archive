# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 17:01:29 2023

@author: bmatt
"""

import numpy as np
import matplotlib.pyplot as plt

n = 1000
##################### CONSTANTS ########################
eta_p = 0.9
T0 = 220
P0 = 80000

gam = 1.4
R = 287
Cp = 1005
hf = 42800000
g = 9.81
pi_c = 18.5

Cd = 0.03195
rho = 1
Sw = 133

Tt4 = 1450

u0_arr = np.linspace(100, 300, n)

fuelWeight = 2500


################### FUNCITONS, do not change!! ########################
def get_M0(u0, gam, T0):
    M0 = u0 / np.sqrt(gam * R * T0)
    return M0

def get_tau_f(hf, Cp, T0):
    tau_f = hf / (Cp* T0)
    return tau_f

def get_tau_c(pi_c, gam):
    tau_c = pi_c **((gam - 1)/gam)
    return tau_c

def get_tau_r(gam, u0, T0):
    tau_r = 1 + (((gam-1)/2)*get_M0(u0, gam, T0)**2)
    return tau_r

def get_tau_gam(Tt4, T0):
    tau_gam = Tt4 / T0
    return tau_gam

def get_f(Tt4, T0, gam, u0, pi_c, hf, Cp):
    num = get_tau_gam(Tt4, T0) - get_tau_r(gam, u0, T0) * get_tau_c(pi_c, gam)
    dom = get_tau_f(hf, Cp, T0) - get_tau_gam(Tt4, T0)
    
    f = num/dom
    return f

def get_tau_t(gam, u0, pi_c, eta_p, Tt4, T0):
    term1 = 1/ (get_tau_r(gam, u0, T0) * get_tau_c(pi_c, gam))
    term2 = (get_tau_r(gam, u0, T0)-1) / (eta_p**2 * get_tau_gam(Tt4, T0))
    
    tau_t = term1+term2
    return tau_t

def get_MoverM(gam, u0, pi_c, eta_p, Tt4, T0):
    num = get_tau_r(gam, u0, T0) * get_tau_c(pi_c, gam) * get_tau_t(gam, u0, pi_c, eta_p, Tt4, T0) - 1
    dom = get_tau_r(gam, u0, T0) - 1
    
    MoverM = np.sqrt(num/dom)
    return MoverM
    
def get_tau_b(gam, u0, pi_c, Tt4, T0):
    tau_b = get_tau_gam(Tt4, T0) / ( get_tau_r(gam, u0, T0) * get_tau_c(pi_c, gam))
    return tau_b

def get_c_core(gam, u0, pi_c, Tt4, T0, eta_p, hf, Cp):
    ToT = np.sqrt(get_tau_b(gam, u0, pi_c, Tt4, T0))
    MoM = get_MoverM(gam, u0, pi_c, eta_p, Tt4, T0)
    f = get_f(Tt4, T0, gam, u0, pi_c, hf, Cp)
    leadingTerm = (gam - 1) * get_M0(u0, gam, T0)**2
    
    c_core = leadingTerm * (((1+f) * MoM * ToT) - 1)
    return c_core

def get_c_fan(Tt4, T0, gam, u0, pi_c, hf, Cp, eta_p):
    f = get_f(Tt4, T0, gam, u0, pi_c, hf, Cp)
    tau_gam = get_tau_gam(Tt4, T0)
    tau_t = get_tau_t(gam, u0, pi_c, eta_p, Tt4, T0)
    tau_r = get_tau_r(gam, u0, T0)
    tau_c = get_tau_c(pi_c, gam)
    
    c_prop = eta_p * ( (1+f)*tau_gam*(1-tau_t) - tau_r*(tau_c - 1) )
    return c_prop


def get_spec_thrust(Tt4, T0, gam, u0, pi_c, hf, Cp, eta_p):
    core_t = get_c_core(gam, u0, pi_c, Tt4, T0, eta_p, hf, Cp) * ((Cp*T0)/(u0))
    prop_t = get_c_fan(Tt4, T0, gam, u0, pi_c, hf, Cp, eta_p) * ((Cp*T0)/(u0))
    
    spec_thrust = core_t + prop_t
    return spec_thrust

def get_drag(u0):
    drag = 0.25 * Cd * rho * u0**2 * Sw
    return drag

def get_range(Tt4, T0, gam, u0, pi_c, hf, Cp, eta_p ):
    spec_t = get_spec_thrust(Tt4, T0, gam, u0, pi_c, hf, Cp, eta_p) 
    f = get_f(Tt4, T0, gam, u0, pi_c, hf, Cp)
    
    rnge = (((fuelWeight) / ((get_drag(u0) / spec_t) * f)) * u0) / 1000
    
    return rnge

##############################################################################

range_arr = np.zeros(n)

for i in range(n):
    range_arr[i] = get_range(Tt4, T0, gam, u0_arr[i], pi_c, hf, Cp, eta_p )
    
plt.plot(u0_arr, range_arr)
plt.grid()
plt.ylabel("Range (km)")
plt.xlabel('Airspeed (m/s)')