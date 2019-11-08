#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 11:06:02 2019

@author: amandaash
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

import numpy as np
import matplotlib.pyplot as plt

g = 9.81

def pendulum_euler(t0, tf, h, l, theta_0, dtheta_0):
    
    C = g/l
    
    theta = theta_0
    dtheta = dtheta_0
    
    x0 = l*np.sin(theta)
    y0 = l*(1-np.cos(theta))
    
    x_array = []
    y_array = []
    
    vx_array = []
    vy_array = []
    
    theta_array = []
    dtheta_array = []
    
    time_array = np.arange(t0,tf,h)
    
    for t in time_array: 
        
        dtheta_n = dtheta - (h*C*np.sin(theta))
        theta_n = theta + (h*dtheta)
        
        x_n = l*np.sin(theta_n)
        y_n = l*(1-np.cos(theta_n))
        
        vx_n = l*np.sin(dtheta_n)
        vy_n = l*(1-np.cos(dtheta_n))
        
        theta_array.append(theta_n)
        dtheta_array.append(dtheta_n)
        x_array.append(x_n)
        y_array.append(y_n)
        vx_array.append(vx_n)
        vy_array.append(vy_n)
        
        theta = theta_n
        dtheta = dtheta_n
        
    return time_array, x_array, y_array, vx_array, vy_array, theta_array, dtheta_array

def pendulum_rk2(t0, tf, h, l, theta_0, dtheta_0):
    
    C = g/l
    
    theta = theta_0
    dtheta = dtheta_0
    
    x0 = l*np.sin(theta)
    y0 = l*(1-np.cos(theta))
    
    x_array = []
    y_array = []
    
    vx_array = []
    vy_array = []
    
    theta_array = []
    dtheta_array = []
    
    time_array = np.arange(t0,tf,h)
    
    for t in time_array:
        #first Euler's method to find half step:
        
        theta_half = theta + (h/2)*dtheta
        dtheta_half = dtheta - (h/2)*C*np.sin(theta)

        dtheta_n = dtheta - h*C*np.sin(theta_half)
        theta_n = theta + (h*dtheta_half)
        
        
        x_n = l*np.sin(theta_n)
        y_n = l*(1-np.cos(theta_n))
        
        vx_n = l*np.sin(dtheta_n)
        vy_n = l*(1-np.cos(dtheta_n))
        
        theta_array.append(theta_n)
        dtheta_array.append(dtheta_n)
        x_array.append(x_n)
        y_array.append(y_n)
        vx_array.append(vx_n)
        vy_array.append(vy_n)
        
        theta = theta_n
        dtheta = dtheta_n
        
    return time_array, vx_array, vy_array, x_array, y_array, theta_array, dtheta_array

def find_period(velocity_array, time_step):
    for index in range(len(list(velocity_array))):
        
        if index == 0 or index == 1:
            continue
        
        elif np.sign(velocity_array[index-1]) != np.sign(velocity_array[index+1]):
            return 2*index*time_step, (2*np.pi)/(2*index*time_step)

degree_radian_conversion = np.pi/180
t0 = 0
tf = 10
h = 0.001
theta_0 = 10*degree_radian_conversion
dtheta_0 = 0
l = 1

t, x, y, vx, vy, angle, d_angle = pendulum_euler(t0, tf, h, l, theta_0, dtheta_0)

period, omega = find_period(d_angle, h)
theory_period = 2*np.pi*np.sqrt(l/g)
print('measured period = {0}s, theoretical period = {1}s'.format(period, theory_period))

plt.title('euler method')
plt.xlabel('time[s]')
plt.ylabel('angle[rad]')
plt.plot(t, angle, '.')
plt.show()

plt.title('euler method')
plt.xlabel('time[s]')
plt.ylabel('$\\frac{d\\theta}{dt}$ [rad $s^{-1}$]')
plt.plot(t, d_angle, '.')
plt.show()

plt.title('euler method')
plt.ylabel('$\\frac{d\\theta}{dt}$ [rad $s^{-1}$]')
plt.xlabel('angle[rad]')
plt.plot(angle, d_angle, '.')
plt.show()


plt.title('euler method')
plt.xlabel('x[m]')
plt.ylabel('y[m]')
plt.plot(x, y, '.')
plt.show()

plt.title('euler method')
plt.xlabel('$V_x$[m $s^{-1}$]')
plt.ylabel('$V_y$[m $s^{-1}$]')
plt.plot(vx, vy, '.')
plt.show()

ax = plt.axes(projection='3d')
ax.set_title('euler method')
ax.set_xlabel('t[s]')
ax.set_ylabel('x[m]')
ax.set_zlabel('y[m]')
ax.scatter3D(t, x, y, color = 'orange', alpha = 0.3)
plt.show()

ax = plt.axes(projection='3d')
ax.set_title('euler method')
ax.set_xlabel('t[s]')
ax.set_ylabel('$V_x$[m $s^{-1}$]')
ax.set_zlabel('$V_y$[m $s^{-1}$]')
ax.scatter3D(t, vx, vy, color = 'purple', alpha = 0.3)
plt.show()

t, x, y,vx, vy, angle, d_angle = pendulum_rk2(t0, tf, h, l, theta_0, dtheta_0)

period, omega = find_period(d_angle, h)
theory_period = 2*np.pi*np.sqrt(l/g)
print('measured period = {0}s, theoretical period = {1}s'.format(period, theory_period))

plt.title('RK2 method')
plt.xlabel('time[s]')
plt.ylabel('angle[rad]')
plt.plot(t, angle, '.')
plt.show()

plt.title('RK2 method')
plt.xlabel('time[s]')
plt.ylabel('$\\frac{d\\theta}{dt}$ [rad $s^{-1}$]')
plt.plot(t, d_angle, '.')
plt.show()

plt.title('RK2 method')
plt.ylabel('$\\frac{d\\theta}{dt}$ [rad $s^{-1}$]')
plt.xlabel('angle[rad]')
plt.plot(angle, d_angle, '.')
plt.show()

plt.title('RK2 method')
plt.xlabel('x[m]')
plt.ylabel('y[m]')
plt.plot(x, y, '.')
plt.show()

plt.title('RK2 method')
plt.xlabel('$V_x$[m $s^{-1}$]')
plt.ylabel('$V_y$[m $s^{-1}$]')
plt.plot(vx, vy, '.')
plt.show()

ax = plt.axes(projection='3d')
ax.set_title('RK2 method')
ax.set_xlabel('t[s]')
ax.set_ylabel('x[m]')
ax.set_zlabel('y[m]')
ax.scatter3D(t, x, y, color = 'orange', alpha = 0.3)
plt.show()

ax = plt.axes(projection='3d')
ax.set_title('RK2 method')
ax.set_xlabel('t[s]')
ax.set_ylabel('$V_x$[m $s^{-1}$]')
ax.set_zlabel('$V_y$[m $s^{-1}$]')
ax.scatter3D(t, vx, vy, color = 'purple', alpha = 0.3)
plt.show()


initial_amp = np.arange(0.001,90,5)
periods = []
for amplitude in initial_amp: 
    theta_0 = amplitude*degree_radian_conversion
    t, x, y,vx, vy, angle, d_angle = pendulum_rk2(t0, tf, h, l, theta_0, dtheta_0)
    period, angular_freq = find_period(d_angle, h)
    periods.append(period)
    
plt.title('RK2 - period')
plt.plot(initial_amp, periods, '.')
plt.axhline(y = 2*np.pi*np.sqrt(l/g))
plt.xlabel('amplitude [m]')
plt.ylabel('period [$s^{-1}$]')
plt.show()
    
        
        
