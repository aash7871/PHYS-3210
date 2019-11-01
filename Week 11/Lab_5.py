#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 19:02:02 2019

@author: amandaash
"""

import numpy as np
import matplotlib.pyplot as plt

def harmonic_oscillator_friction_beta(p,k,v0,x0,m,time_step,t0,tf, mu_s, mu_k, b):
    g = 9.81
    v = v0
    x = x0
    x_val = []
    v_val = []
    time_array = np.arange(t0,tf, time_step)
    t_val = []
    
    
    for n in time_array:
        
        if np.isclose(v,0, atol = 10**-5):
       
            a = (-k/m)*x**(p-1) + (mu_s*g)
        
            if (mu_s*m*g) >= np.abs((k)*x**(p-1)):
                print('Static Friction failure, F_f = {0} N, F_spring = {1} N'.format(str((mu_s*m*g)), str(np.abs((k)*x**(p-1)))))
                break
            else:
                t_val.append(n)
                x_half = x + (time_step/2)*v
                v_half = v + ((time_step/2)*(a))
                a_half = (-k/m)*x_half**(p-1) + (mu_s*g)
        
                vf = v + ((time_step)/m)*(a_half)
                xf = x + (time_step)*v_half
        
                x_val.append(x)
                v_val.append(v)
        
        
                v = vf
                x = xf
            
        else:
            #print('False')
            a = (-k/m)*x**(p-1) - (mu_k*g*(v/np.abs(v))) - ((b/m)*v)
            t_val.append(n)
        #first Euler's method to find half step:
            x_half = x + (time_step/2)*v
            v_half = v + ((time_step/2)*(a))
            a_half = (-k/m)*x_half**(p-1) - (mu_k*g*(v_half/np.abs(v_half))) - ((b/m)*v_half)
        
            vf = v + ((time_step)/m)*(a_half)
            xf = x + (time_step)*v_half
        
            x_val.append(x)
            v_val.append(v)
        
        
            v = vf
            x = xf
    
    return x_val, v_val, t_val

def harmonic_oscillator_drive(p,k,v0,x0,m,time_step,t0,tf, F0, omega, mu_s = 0, mu_k = 0, b = 0):
    v = v0
    x = x0
    x_val = []
    v_val = []
    time_array = np.arange(t0,tf, time_step)

    for n in time_array:
        a = ((-k/m)*x**(p-1))+((F0/m)*np.sin(omega*n)) 
        #first Euler's method to find half step:
        x_half = x + (time_step/2)*v
        v_half = v + ((time_step/2))*a
        a_half = ((-k/m)*x_half**(p-1))+((F0/m)*np.sin(omega*(n+(time_step/2))))
        
        vf = v + ((time_step)/m)*(a_half)
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
            return 2*index*time_step, (2*np.pi)/(2*index*time_step)

dt = 0.0001
mass = 1
p_value = 2
k_constant = 100
v_initial = 0
x_initial = 1
t_initial = 0
t_final = 10
static_coeff = 0.45
kinetic_coeff = 0.35
viscous_coeff = 0.6

plt.title('damped oscillator, P = {0}, k = {1}, $\\mu_s$ = {2}, $\\mu_k$ = {3}, b = {4}' .format(p_value, k_constant, static_coeff, kinetic_coeff, viscous_coeff))
x_val,v_val,t_val = harmonic_oscillator_friction_beta(p_value,k_constant,v_initial,x_initial,mass,dt,t_initial,t_final,static_coeff,kinetic_coeff,viscous_coeff)
period, angular_frequency = find_period(v_val, dt)
#print(angular_frequency, angular_frequency*2*m)
plt.plot(x_val, t_val)
plt.xlabel('x[m]')
plt.ylabel('t[s]')
#plt.plot(v_val, t_val)
plt.show()


dt = 0.0001
mass = 1
p_value = 2
k_constant = 1
v_initial = 0
x_initial = 1
t_initial = 0
t_final = 100
F_drive = 10000
frequency_drive = 10

#Large Driving Force: 

plt.title('driven oscillator, P = {0}, k = {1}, $F_0$ = {2}, $\\omega$ = {3}'.format(p_value, k_constant, F_drive, frequency_drive))
x_drive, v_drive, t_drive = harmonic_oscillator_drive(p_value,k_constant,v_initial,x_initial,mass,dt,t_initial,t_final,F_drive,frequency_drive)
plt.plot(x_drive, t_drive, '-')
plt.xlabel('x[m]')
plt.ylabel('t[s]')
plt.show()

#beats conditions?: dt = 0.0001, m = 1, p = 2, k = 10, v0 = 0, x0 = 1, t0 = 0, tf = 10, F0 = 10, omega = 1

dt = 0.0001
mass = 1
p_value = 4
k_constant = 10
v_initial = 0
x_initial = 1
t_initial = 0
t_final = 10
F_drive = 10
frequency_drive = 1

plt.title('driven oscillator, P = {0}, k = {1}, $F_0$ = {2}, $\\omega$ = {3}'.format(p_value, k_constant, F_drive, frequency_drive))
x_drive, v_drive, t_drive = harmonic_oscillator_drive(p_value,k_constant,v_initial,x_initial,mass,dt,t_initial,t_final,F_drive,frequency_drive)
plt.plot(x_drive, t_drive, '-')
plt.xlabel('x[m]')
plt.ylabel('t[s]')
plt.show()


#resonance conditions?: dt = 0.001, m = 1, p = 2, k = 1, v0 = 0, x0 = 1, t0 = 0, tf = 40, F0 = 1, omega = 1


