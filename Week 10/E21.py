#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 09:59:35 2019

@author: amandaash
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

def harmonic_oscillator(p,k,v0,x0,m,time_step,t0,tf):
    v = v0
    x = x0
    x_val = []
    v_val = []
    time_array = np.arange(t0,tf, time_step)

    for n in time_array:
        #first Euler's method to find half step:
        x_half = x + (time_step/2)*v
        v_half = v + ((time_step/2)*(-k/m)*x**(p-1))
        
        vf = v + ((time_step)/m)*(-k*x_half**(p-1))
        xf = x + (time_step)*v_half
        
        x_val.append(x)
        v_val.append(v)
        
        v = vf
        x = xf
    
    return x_val, v_val, time_array

def find_period(velocity_array, time_step):
    for index in range(len(list(velocity_array))):
        
        if index == 0 or index == 1:
            continue
        
        elif np.sign(velocity_array[index-1]) != np.sign(velocity_array[index+1]):
            return 2*index*time_step


dt = 0.00001
x_values, v_values, times = harmonic_oscillator(6, 10, 0, 1, 1, dt, 0, 10)     
period_index = find_period(v_values, dt)
plt.plot(times, x_values, '-')
plt.axvline(x = period_index)
plt.show()
print(period_index)

periods = []
amplitudes = np.arange(1, 5, 0.25)
for amp in amplitudes:
    x_values, v_values, times = harmonic_oscillator(6, 10, 0, amp, 1, dt, 0, 10)
    period = find_period(v_values, dt)
    periods.append(period)
    """
    plt.plot(x_values, times, '-')
    plt.axhline(y = period, color = 'k')

    plt.xlabel('x[m]')
    plt.ylabel('t[s]')
    plt.show()
    """

plt.plot(amplitudes, periods, '.')
plt.xlabel('A[m]')
plt.ylabel('T[s]')
plt.show()


P_array = np.arange(2,14,2)
period_nonharmonic = []
k = 10
m = 1
T_theory = 2*np.pi*np.sqrt(m/k)
for P in P_array:
    x,v,t = harmonic_oscillator(P, k, 0, 1, m, dt, 0, 10)
    period = find_period(v, dt)
    oscillations50 = period*50
    #changed my timestep because computation time
    x50,v50,t50 = harmonic_oscillator(P, k, 0, 1, m, 0.01, 0, oscillations50)
    x50 = np.array(x50)
    v50 = np.array(v50)
    t50 = np.array(t50)
    PE = 0.5*k*x50**2
    KE = 0.5*m*v50**2
    E_tot = PE+KE
    
    print(np.mean(KE), (P/2)*np.mean(PE))
    
    fig1 = plt.figure(figsize = (25,5))
    plt.title('P = {0}'.format(P))
    plt.plot(t50, PE, '-', label = 'potential energy')
    plt.plot(t50, KE, '-', label = 'kinetic energy')
    plt.plot(t50, E_tot, '-', label = 'total mechanical energy')
    plt.legend()
    plt.xlabel('t[s]')
    plt.ylabel('E[J]')
    plt.savefig('P_{0}.pdf'.format(P))
    fig1.show()
    
    fig2 = plt.figure()
    stability = -np.log10(np.abs((E_tot-E_tot[0])/E_tot[0]))
    plt.plot(t50, stability, '.')
    plt.xlabel('t[s]')
    plt.ylabel('$log(\\frac{E_{tot}-E_{t = 0}}{E_{t=0}})$')
    plt.savefig('error_P{0}.pdf'.format(P))
    plt.show()
    
    period_nonharmonic.append(period)

plt.plot(P_array, period_nonharmonic, '.')
plt.xlabel('P')
plt.ylabel('T[s]')
plt.axhline(y = T_theory)
plt.show()


    

    
    