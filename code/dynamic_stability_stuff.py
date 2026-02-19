# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 19:38:45 2023

@author: bmatt
"""

import numpy as np
import matplotlib.pyplot as plt


n = 1000
g = 9.81
S = 27.5
rho = 1
CD_u = -0.03 
CL_0 = 0.386    
CD_0 = 0.03105
m = 7000
u = 270
gam = 1.4
R = 287
T = 220

########### CURRENTLY NEED THESE NUMBERS ###########
CL_a = 2.34
Cmq = -3.037
cbar = 2.5
Iy = -751289.3809
Cma = 1.496
Cma_dot = -1.382
#####################################################

def get_CLu(u, CL_0):
    M = u / np.sqrt(gam*R*T)
    CLu = -((M**2) / (1-M**2)) * CL_0
    return CLu

def get_CDu(u, CD_0):
    M = u / np.sqrt(gam*R*T)
    CDu = -((M**2) / (1-M**2)) * CD_0
    return CDu

def get_Q(rho, u):
    Q = 0.5*rho*u**2
    return Q

def get_Zu(rho, u, CL_0, S, m):
    num = -(get_CLu(u, CL_0) + 2*CL_0) * S * get_Q(rho, u)
    dom = m*u
    Zu = num/dom
    return Zu

def get_Xu(CD_0, S, m, u, rho):
    num = -(get_CDu(u, CD_0) + 2* CD_0) * S * get_Q(rho, u)
    dom = m*u
    assert dom != 0, 'divide by zero'
    Xu = num / dom
    return Xu

def get_phugoid_freq(rho, u, CL_0, S, m):
    num= get_Zu(rho, u, CL_0, S, m)*g
    dom = u
    assert dom != 0, 'divide by zero, phug freq'
    assert num/dom >= 0, 'imaginary root, phug freq'
    phugoid_freq = np.sqrt(num/dom)
    return phugoid_freq

def get_phugoid_damp(rho, u, CL_0, S, m, CD_0):
    num = -get_Xu(CD_0, S, m, u, rho)
    dom = 2* get_phugoid_freq(rho, u, CL_0, S, m)
    assert dom != 0, 'divide by zero, phug damp'
    phugoid_damp = num / dom
    return phugoid_damp


print('phugoid freq = ' + str(get_phugoid_freq(rho, u, CL_0, S, m)))
print('phugoid damp = ' + str(get_phugoid_damp(rho, u, CL_0, S, m, CD_0)))


def get_Zw(CL_a, CD_0, rho, u, S, m):
    num = -(CL_a + CD_0) * S * get_Q(rho, u)
    dom = m*u
    Zw = num / dom
    assert dom != 0, 'divide by zero in Zw'
    return Zw

def get_Za(CL_a, CD_0, rho, u, S, m):
    Za = u * get_Zw(CL_a, CD_0, rho, u, S, m)
    return Za

def get_Mq(Cmq, cbar, u, S, rho, Iy):
    num = Cmq*cbar**2*S*get_Q(rho, u)
    dom = 2*u*Iy
    Mq = num/dom
    assert dom != 0, 'divide by zero, Mq'
    return Mq

def get_Mw(Cma, S, rho, u, cbar, Iy):
    num = Cma*S*cbar*get_Q(rho, u)
    dom = u * Iy
    Mw = num / dom
    assert dom != 0, 'divide by zero, Mw'
    return Mw

def get_Ma(Cma, S, rho, u, cbar, Iy):
    Ma = get_Mw(Cma, S, rho, u, cbar, Iy) * u
    return Ma

def get_Mw_dot(Cma_dot, cbar, rho, u, S, Iy):
    num = Cma_dot * cbar * get_Q(rho, u) * S * cbar
    dom = 2 * u * u * Iy
    Mw_dot = num / dom
    return Mw_dot

def get_Ma_dot(Cma_dot, cbar, rho, u, S, Iy):
    Ma_dot = u * get_Mw_dot(Cma_dot, cbar, rho, u, S, Iy)
    return Ma_dot


def get_sp_freq(CL_a, CD_0, rho, u, S, m, Cmq, cbar, Iy, Cma):
    num = get_Za(CL_a, CD_0, rho, u, S, m)*get_Mq(Cmq, cbar, u, S, rho, Iy)
    dom = u
    term1 = num / dom
    sp_freq = np.sqrt(term1 - get_Ma(Cma, S, rho, u, cbar, Iy))
    return sp_freq

def get_sp_damp(Cmq, cbar, u, S, rho, Iy, Cma_dot, Cma, CL_a):
    Mq = get_Mq(Cmq, cbar, u, S, rho, Iy)
    Ma_dot = get_Ma_dot(Cma_dot, cbar, rho, u, S, Iy)
    Za = get_Za(CL_a, CD_0, rho, u, S, m)
    freq = get_sp_freq(CL_a, CD_0, rho, u, S, m, Cmq, cbar, Iy, Cma)
    
    num = Mq + Ma_dot + (Za / u)
    dom = 2 * freq
    sp_damp = - (num / dom)
    
    return sp_damp

print('sp freq = ' +str(get_sp_freq(CL_a, CD_0, rho, u, S, m, Cmq, cbar, Iy, Cma)))
print('sp damp = ' +str(get_sp_damp(Cmq, cbar, u, S, rho, Iy, Cma_dot, Cma, CL_a)))                             



    

