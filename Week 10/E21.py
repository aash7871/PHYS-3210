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

def sine_fit(time, A, phi, omega):
    y = A*np.sin((omega*time)+phi)
    return y

amplitudes = np.arange(1, 1.5, 0.25)
for amp in amplitudes:
    x_values, v_values, times = harmonic_oscillator(6, 10, 0, amp, 1, 0.0001, 0, 10)
    parameters, parameters_covariance =  opt.curve_fit(sine_fit, times, x_values)
    print(parameters)
    plt.plot(x_values, times)
    plt.plot(sine_fit(parameters[0], parameters[1], parameters[2], times), times)
    plt.show()