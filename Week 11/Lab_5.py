#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 19:02:02 2019

@author: amandaash
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from matplotlib.animation import FuncAnimation

def harmonic_oscillator_friction(p,k,v0,x0,m,time_step,t0,tf, mu_s, mu_k, b):
    g = 9.81
    v = v0
    x = x0
    x_val = []
    v_val = []
    time_array = np.arange(t0,tf, time_step)
    t_val = []
    
    
    for n in time_array:
        
        if np.isclose(v,0, atol = time_step):
       
            a = (-k/m)*x**(p-1) + (mu_s*g)
        
            if (mu_s*m*g) >= np.abs((k)*x**(p-1)):
                print('Static Friction failure, F_f = {0}, F_spring = {1}'.format(str((mu_s*m*g)), str(np.abs((k)*x**(p-1)))))
                break
            else:
                t_val.append(n)
            
        else:
            #print('False')
            a = (-k/m)*x**(p-1) - (mu_k*g*(v/np.abs(v))) - ((b/m)*v)
            t_val.append(n)
        #first Euler's method to find half step:
        x_half = x + (time_step/2)*v
        v_half = v + ((time_step/2)*(a))
        
        vf = v + ((time_step)/m)*(a)
        xf = x + (time_step)*v_half
        
        x_val.append(x)
        v_val.append(v)
        
        
        v = vf
        x = xf
    
    return x_val, v_val, t_val

def find_period(velocity_array, time_step):
    for index in range(len(list(velocity_array))):
        
        if index == 0 or index == 1:
            continue
        
        elif np.sign(velocity_array[index-1]) != np.sign(velocity_array[index+1]):
            return 2*index*time_step, (2*np.pi)/(2*index*time_step)

dt = 0.001
m = 1
x_val,v_val,t_val = harmonic_oscillator_friction(2,10,0,1,1,dt,0,10,0.15,0.1,0)
period, angular_frequency = find_period(v_val, dt)
#print(angular_frequency, angular_frequency*2*m)
plt.plot(x_val, t_val)
#plt.plot(v_val, t_val)
plt.show()

