#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 16:59:06 2019

@author: amandaash
"""

import numpy as np
import matplotlib.pyplot as plt

def E_laser_euler(t0, tf, h, E0, dE0, omega, tau, g, g_tild):
    time_array = np.arange(t0, tf, h)
    E = E0
    dE = dE0
    E_array = []
    dE_array = []
    for t in time_array:
        
        alpha = (-(omega**2)*E)-(dE/tau)+((g-(g_tild*E**2))*dE)
        
        E = E + (h*dE)
        
        dE = dE + (h*alpha)
        
        E_array.append(E)
        
        dE_array.append(dE)
        
    return E_array, dE_array, time_array

def E_laser_RK2(t0, tf, h, E0, dE0, omega, tau, g, g_tild):
    time_array = np.arange(t0, tf, h)
    E = E0
    dE = dE0
    E_array = []
    dE_array = []
    for t in time_array:
        
        alpha = (-(omega**2)*E)-(dE/tau)+((g-(g_tild*E**2))*dE)
        
        E_half = E + ((h/2)*dE)
        
        dE_half = dE + ((h/2)*alpha)
        
        alpha_half = (-(omega**2)*E_half)-(dE_half/tau)+((g-(g_tild*E_half**2))*dE_half)
        
        E = E + ((h)*dE_half)
        
        dE = dE + ((h)*alpha_half)
        
        E_array.append(E)
        
        dE_array.append(dE)
        
    return E_array, dE_array, time_array

t0 = 0
tf = 50
h = 1e-4
omega = 1
tau = -1000
g = 0.1
g_tild = 0
E0 = 5
dE0 = 0

E_val, dE_val, time = E_laser_euler(t0, tf, h, E0, dE0, omega, tau, g, g_tild)
plt.title('E')
plt.plot(time, E_val)
plt.show()

"""
plt.title('dE')
plt.plot(time, dE_val)
plt.show()
"""

plt.plot(E_val, dE_val)
plt.show()

t0 = 0
tf = 100
h = 1e-4
omega = 0.5
tau = 1
g = 1
g_tild = 0.1
E0 = 0.5
dE0 = 0

E_val, dE_val, time = E_laser_euler(t0, tf, h, E0, dE0, omega, tau, g, g_tild)

plt.title('E')
plt.plot(time, E_val)
plt.show()
"""
plt.title('dE')
plt.plot(time, dE_val)
plt.show()
"""
plt.plot(E_val, dE_val)
plt.show()


################################

t0 = 0
tf = 200
h = 1e-4
omega = 0.5
tau = 1
g = 1
g_tild = 0
E0 = 1
dE0 = 0

t0 = 0
tf = 50
h = 1e-4
omega = 1
tau = 1.5
g = 1
g_tild = 0
E0 = 5
dE0 = 0

E_val, dE_val, time = E_laser_RK2(t0, tf, h, E0, dE0, omega, tau, g, g_tild)
plt.title('E')
plt.plot(time, E_val)
plt.xlabel('t')
plt.ylabel('E(t)')
plt.tight_layout()
plt.savefig('g_tilde_0.pdf')
plt.show()

"""
plt.title('dE')
plt.plot(time, dE_val)
plt.show()
"""

plt.plot(E_val, dE_val)
plt.show()

t0 = 0
tf = 100
h = 1e-4
omega = 0.5
tau = 2
g = 1
g_tild = 0.5
E0 = 0.5
dE0 = 0

E_val, dE_val, time = E_laser_RK2(t0, tf, h, E0, dE0, omega, tau, g, g_tild)

plt.plot(time, E_val)
plt.xlabel('t')
plt.ylabel('E(t)')
plt.tight_layout()
plt.savefig('g_tilde_nonzero.pdf')
plt.show()
"""
plt.title('dE')
plt.plot(time, dE_val)
plt.show()
"""
plt.plot(E_val, dE_val)
plt.xlabel('E')
plt.ylabel("$\dot{E}$")
plt.tight_layout()
plt.savefig('phase_space_tilde_nonzero.pdf')
plt.show()