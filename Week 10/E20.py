#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 09:56:21 2019

@author: amandaash
"""

import numpy as np
import matplotlib.pyplot as plt


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
        #v = v + ((n + (time_step/2))/m)*(-k*x1**(p-1))
        v = v + ((time_step)/m)*(-k*x_half**(p-1))
        x = x + (time_step)*v_half
        
        x_val.append(x)
        v_val.append(v)
    
    return x_val, v_val, time_array

P_val = np.arange(2,14,2)
fig1 = plt.figure()
ax1 = fig1.add_subplot(1,1,1)
fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)
for P in P_val:

    x,v,t = harmonic_oscillator(P,1,1,0,1,0.0001,0,10)
    ax1.plot(t,x, label = 'P = {0}'.format(P))
    ax2.plot(t,v, label = 'P = {0}'.format(P))

ax1.legend()
fig1.show()


    