#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 22:35:20 2019

@author: amandaash
"""

import numpy as np
import matplotlib.pyplot as plt

def acceleration():

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
        vx_n = (((h*y**2)/m)*((2*x*np.exp(-(x**2+y**2)))-(2*(x**3)*np.exp(-(x**2+y**2))))) + vx
        vy_n = (((h*x**2)/m)*((2*y*np.exp(-(x**2+y**2)))-(2*(y**3)*np.exp(-(x**2+y**2))))) + vy
        
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

t, x_t, y_t, vx_t, vy_t = scatter(0,0,1,1,0,10,0.001,10)

plt.plot(x_t,y_t)
plt.show()  