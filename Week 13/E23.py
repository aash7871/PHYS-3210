#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 22:35:20 2019

@author: amandaash
"""

import numpy as np
import matplotlib.pyplot as plt

def acceleration(x, y, m):
    ax = -((2*x*y**2)/m)*np.exp((-x**2)-(y**2))*(1-x**2)
    ay = -((2*y*x**2)/m)*np.exp((-x**2)-(y**2))*(1-y**2)
    
    return ax, ay

def scatter(x0,y0,vx0,vy0,t0,tf,h,m):
    x = x0
    y = y0
    vx = vx0
    vy = vy0
    
    
    time_array = np.arange(t0,tf,h)
    
    vx_val = []
    vy_val = []
    
    x_val = []
    y_val = []
     
    for step in time_array:

        if np.abs(x)>=3:
            vx = -vx
            vy = vy
            x = x
            y = y
            vx_val.append(vx)
            vy_val.append(vy)
        
            x_val.append(x)
            y_val.append(y)
    
        if np.abs(y)>=3:
            vx = vx
            vy = -vy
            x = x
            y = y
            vx_val.append(vx)
            vy_val.append(vy)
        
            x_val.append(x)
            y_val.append(y)
        
        else:
            ax, ay = acceleration(x,y,m)
        
        #vx_n = (((h*y**2)/m)*((2*x*np.exp(-(x**2+y**2)))-(2*(x**3)*np.exp(-(x**2+y**2))))) + vx
        #vy_n = (((h*x**2)/m)*((2*y*np.exp(-(x**2+y**2)))-(2*(y**3)*np.exp(-(x**2+y**2))))) + vy
            vx_n = h*ax + vx
            vy_n = h*ay + vy
        
            x_n = (h*vx) + x
            y_n = (h*vy) + y
        
            vx_val.append(vx_n)
            vy_val.append(vy_n)
        
            x_val.append(x_n)
            y_val.append(y_n)
        
            x = x_n
            y = y_n
        
            vx = vx_n
            vy = vy_n
        
    return time_array, x_val, y_val, vx_val, vy_val

def potential(x,y):
    potential_function = x**2*y**2*np.exp(-x**2-y**2)
    return potential_function

x_grid = np.arange(-15, 15, 0.01)
y_grid = np.arange(-15, 15, 0.01)
X_mesh, Y_mesh = np.meshgrid(x_grid, y_grid)
V = potential(X_mesh, Y_mesh)

fig, ax = plt.subplots(1, 1, figsize=(10, 8))
t1, x_t1, y_t1, vx_t1, vy_t1 = scatter(0.25,0.5,0,-0.1,0,150,0.001,0.1)
t2, x_t2, y_t2, vx_t2, vy_t2 = scatter(-0.25,0.5,0,-0.1,0,150,0.001,0.1)


cb = ax.contourf(X_mesh, Y_mesh, V, 10)
ax.set_ylim(-3,3)
ax.set_xlim(-3,3)
ax.plot(x_t1, y_t1)
ax.plot(x_t2, y_t2)
plt.show()


t_wild, x_wild, y_wild, vx_wild, vy_wild = scatter(-0.2,0.9,0,-0.01,0,1000,0.01,0.1)
plt.plot(x_wild,y_wild)
#plt.ylim(-3,3)
#plt.xlim(-3,3)
plt.show() 

 

