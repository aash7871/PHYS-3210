#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 09:57:43 2019

@author: amandaash
"""

import numpy as np
import matplotlib.pyplot as plt

p = 2
v = 1
x = 0
m = 10
time_step = 0.0001
k = 3
t0 = 0
tf = 10

"""
x_val = []
v_val = []

time_array = np.arange(t0,tf, time_step)

for n in time_array:
    
    v1 = v + (time_step/m)*(-k*x**(p-1))
    x1 = x + time_step*v
    
    x_val.append(x1)
    v_val.append(v1)
    
    v1 = v
    x1 = x

plt.plot(time_array, x_val) 
plt.show()
plt.plot(time_array, v_val)
plt.show()
"""

def harmonic_oscillator(p,k,v0,x0,m,time_step,t0,tf):
    v = v0
    x = x0
    x_val = []
    v_val = []
    time_array = np.arange(t0,tf, time_step)

    for n in time_array:
    
        vf = v + (time_step/m)*(-k*x**(p-1))
        xf = x + time_step*v
    
        x_val.append(xf)
        v_val.append(vf)
        
        x = xf
        v = vf
    
    return x_val, v_val, time_array


#P_val = np.arange(2,8,2)
P_val = np.array([2,6,10])
fig1 = plt.figure()  
ax1 = fig1.add_subplot(1, 1, 1)
fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)

for P_value in P_val:
       x_P, v_P, t_P= harmonic_oscillator(P_value, 10, 0, 1, 1, 0.0001, 0, 10)
       ax1.plot(x_P, t_P, label = "P = {0}".format(P_value))
       ax2.plot(v_P, t_P, label = "P = {0}".format(P_value))

ax1.set_xlabel('distance')   
ax1.set_ylabel('time')       
ax1.legend()
fig1.savefig('spring_dt.pdf')
fig1.show()

ax2.set_xlabel('velocity')   
ax2.set_ylabel('time')   
ax2.legend()
fig2.savefig('spring_vt.pdf')
fig2.show()

#amplitude - frequency things: 
fig3 = plt.figure()
ax3 = fig3.add_subplot(1,1,1)
x_ic = np.arange(0.5,2.0,0.5)
for amplitude in x_ic: 
    x_a, v_a, t_a = harmonic_oscillator(6, 10, 0, amplitude, 1, 0.0001, 0, 10)
    ax3.plot(x_a, t_a, label = '$x_0$ = {0}'.format(amplitude))
ax3.set_title('P = 6, non-harmonic oscillator varying $x_0$')
ax3.set_xlabel('x')
ax3.set_ylabel('t')
ax3.legend()
fig3.savefig('non_harmonic_amplitude.pdf')
fig3.show()


#Going between the RK2 method and the Euler method from exercise 19, I see no 
#difference between the two methods in either the position vs. time or velocity vs. time for the oscillator. 