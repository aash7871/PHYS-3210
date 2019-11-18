#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 22:35:20 2019

@author: amandaash
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

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
            x = x + (vx*h)
            y = y +(vy*h)
            vx_val.append(vx)
            vy_val.append(vy)
        
            x_val.append(x)
            y_val.append(y)
    
        if np.abs(y)>=3:
            vx = vx
            vy = -vy
            x = x+(vx*h)
            y = y+(vy*h)
            vx_val.append(vx)
            vy_val.append(vy)
        
            x_val.append(x)
            y_val.append(y)
        
        else:
            ax, ay = acceleration(x,y,m)
    
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

def scatter_RK2(x0,y0,vx0,vy0,t0,tf,h,m):
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
            x = x + (vx*h)
            y = y +(vy*h)
            vx_val.append(vx)
            vy_val.append(vy)
        
            x_val.append(x)
            y_val.append(y)
    
        if np.abs(y)>=3:
            vx = vx
            vy = -vy
            x = x+(vx*h)
            y = y+(vy*h)
            vx_val.append(vx)
            vy_val.append(vy)
        
            x_val.append(x)
            y_val.append(y)
        
        else:
            
            x_half = x + (time_step/2)*v
            v_half = v + ((time_step/2)*(-k/m)*x**(p-1))
        
            vf = v + ((time_step)/m)*(-k*x_half**(p-1))
            xf = x + (time_step)*v_half
        
            x_val.append(x)
            v_val.append(v)
            
            ax, ay = acceleration(x,y,m)
            
            x_half = x + (h/2)*vx
            y_half = y + (h/2)*vy
            
            vx_half = vx + ((h/2)*(ax))
            vy_half = vy + ((h/2)*(ay))
    
            vx_n = h*ax + vx_half
            vy_n = h*ay + vy_half
        
            x_n = (h*vx_half) + x_half
            y_n = (h*vy_half) + y_half
        
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
t1, x_t1, y_t1, vx_t1, vy_t1 = scatter(0.25,0.5,0,-0.1,0,200,0.001,0.1)
t2, x_t2, y_t2, vx_t2, vy_t2 = scatter(-0.25,0.5,0,-0.1,0,200,0.001,0.1)


cb = ax.contourf(X_mesh, Y_mesh, V, 10)
ax.set_ylim(-3,3)
ax.set_xlim(-3,3)
ax.plot(x_t1, y_t1)
ax.plot(x_t2, y_t2)
plt.show()

fig, ax = plt.subplots(1, 1, figsize=(10, 8))
t_wild, x_wild, y_wild, vx_wild, vy_wild = scatter(-0.2,0.94999,0,-1,0,100,0.01,0.1)
cb = ax.contourf(X_mesh, Y_mesh, V, 10)
ax.set_ylim(-3,3)
ax.set_xlim(-3,3)
ax.plot(x_wild, y_wild)
plt.show()


x_grid = np.arange(-3, 3, 0.01)
y_grid = np.arange(-3, 3, 0.01)
X_mesh, Y_mesh = np.meshgrid(x_grid, y_grid)
V = potential(X_mesh, Y_mesh)
fig = plt.figure(figsize = (10,10))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(-3,3)
ax.set_ylim(-3,3)
ax.plot3D(x_t1, y_t1, potential(np.array(x_t1), np.array(y_t1)))
ax.plot3D(x_t2, y_t2, potential(np.array(x_t2), np.array(y_t2)))
ax.plot_wireframe(X_mesh, Y_mesh, V, color='purple', alpha = 0.25,rstride=5, cstride=5)

#plt.show()

 

