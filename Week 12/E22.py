#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 09:58:52 2019

@author: amandaash
"""

import numpy as np
import matplotlib.pyplot as plt
#x_match should be some very large number

"""
def wave_funcn(psi_0, dpsi_0, x0, xf, h,  m, x_match):
    psi = psi_0
    dpsi = dpsi_0
    hbar = (6.626*10**-34)/(2*np.pi)
    k = 2*m/hbar**2
    
    x_array = np.arange(x0,xf,h)
    
    psi_val = []
    dpsi_val = []
    
    for x in x_array:
        # potential energy conditions
        if abs(x) <= 2:

            V = -83
            
            psi_n = h*dpsi + psi
            dpsi_n = h*((psi*k**2)+ (((2*m)/hbar**2)*V*psi)) + dpsi
        
            psi_val.append(psi_n)
            dpsi_val.append(dpsi_n)
        
            psi = psi_n
            dspi = dpsi_n
        
        else:
            if x < 2
            V = 0
            
            psi_n = 0
            dpsi_n = 0
            
            psi_val.append(psi_n)
            dpsi_val.append(dpsi_n)
        
            psi = psi_n
            dspi = dpsi_n
            
        
        
    return x_array, psi_val, dpsi_val
    
"""


def wave_funcn_alpha(x0, xf, h, a, E, m, V_well):

    x_array = np.arange(x0, xf, h)
    kappa = -C*E
    
    x_well = []
    psi_val_out = []
    dpsi_val_out = []
    
    x_out = []
    psi_val_well = []
    dpsi_val_well = []
    
    
    for x in x_array: 
        if abs(x) >= a: 
            x_out.append(x)
            if x < 0: 
                psi = np.exp(x*kappa)
                dpsi = (kappa)* psi
                
                
                psi_val_out.append(psi)
                dpsi_val_out.append(dpsi)
            if x > 0:
                
                psi = np.exp(-x*kappa)
                dpsi = -(kappa)* psi
                
                psi_val_out.append(psi)
                dpsi_val_out.append(dpsi)
        
        else: 
            x_well.append(x)
            psi_n = h*dpsi + psi
            dpsi_n = h*((psi*kappa)+ ((C*V_well*psi))) + dpsi
        
            psi_val_well.append(psi_n)
            dpsi_val_well.append(dpsi_n)
        
            psi = psi_n
            dpsi = dpsi_n
        
    return x_well, x_out, psi_val_well, dpsi_val_well, psi_val_out, dpsi_val_out, 

x0 = -5
xf = 5
h = 0.001
a = 2
E = -16
#m = 1.67*10**-27
m = 9.11*10**-31
C = 0.0483
V = -83


x_well_array, x_out_array, psi_array_well, dpsi_array_well, psi_array_out, dpsi_array_out = wave_funcn_alpha(x0,xf,h,a,E,m, V)
plt.plot(x_well_array, psi_array_well, '.')
plt.plot(x_out_array, psi_array_out, '.')
plt.show()
        
    
    