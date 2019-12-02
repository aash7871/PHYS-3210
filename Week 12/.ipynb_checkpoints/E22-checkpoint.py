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
    
    psi_val_out_l = []
    dpsi_val_out_l = []
    
    psi_val_out_r = []
    dpsi_val_out_r = []
    
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
                
                psi_val_out_l.append(psi)
                dpsi_val_out_l.append(dpsi)
            if x > 0:
                
                psi = np.exp(-x*kappa)
                dpsi = -(kappa)* psi
                
                psi_val_out_r.append(psi)
                dpsi_val_out_r.append(dpsi)
                
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
    
    
    if x0 > xf: 
        #coming from the right
        #print(dpsi_val_out_l[0], psi_val_out_l[0])
        zl = dpsi_val_well[-1]/psi_val_well[-1]
        zr = dpsi_val_out_l[0]/psi_val_out_l[0]
        z = (zl - zr)/(zl+zr)
    if x0 < xf: 
        #coming from the left
        #print(dpsi_val_out_r[0], psi_val_out_r[0])
        zr = dpsi_val_well[-1]/psi_val_well[-1]
        zl = dpsi_val_out_r[0]/psi_val_out_r[0]
        z = (zl - zr)/(zl+zr)
        
       
    return x_well, x_out, psi_val_well, dpsi_val_well, psi_val_out, dpsi_val_out, z

x0 = -3
xf = 3
h = 0.0001
a = 2
E = -10
m = 1.67*10**-27
C = 0.0483
V = -83


x_well_array, x_out_array, psi_array_well, dpsi_array_well, psi_array_out, dpsi_array_out, z = wave_funcn_alpha(x0,xf,h,a,E,m, V)
plt.plot(x_well_array, psi_array_well, '.', ms = 0.5)
plt.plot(x_out_array, psi_array_out, '.', ms = 0.5)
plt.show()
        
#energy solver 
#going from left and right for a single energy:
E0 = -10
E1 = -20

for n in range(20):
    x_well_array, x_out_array, psi_array_well, dpsi_array_well, psi_array_out, dpsi_array_out, z0 = wave_funcn_alpha(x0,xf,h,a,E0,m, V) 
    print(z0)
    #x_well_array, x_out_array, psi_array_well, dpsi_array_well, psi_array_out, dpsi_array_out, z0l = wave_funcn_alpha(xf,x0,-h,a,E0,m, V) 
    #print(z0l, z0r)
    #delta_E0 = (z0l - z0r)/(z0l + z0r)
    #print(delta_E0)
    
    x_well_array, x_out_array, psi_array_well, dpsi_array_well, psi_array_out, dpsi_array_out, z1 = wave_funcn_alpha(x0,xf,h,a,E1,m, V) 
    print(z1)
    #x_well_array, x_out_array, psi_array_well, dpsi_array_well, psi_array_out, dpsi_array_out, z1l = wave_funcn_alpha(xf,x0,-h,a,E1,m, V)

    #print(z1l, z1r)
    #delta_E1 = (z1l - z1r)/(z1l + z1r)
    #print(delta_E1)
    
    m = (z1 - z0)/(E1-E0)
    En = E0 - (z0/m)
    print(En)
    E0 = En
    E1 = En + (np.sign(z0*m)*7)
  
    #print(En, delta_E0)
    if np.abs(z0*m) <= 1e-4:
        E_match = En
        x_well_array, x_out_array, psi_array_well, dpsi_array_well, psi_array_out, dpsi_array_out, z = wave_funcn_alpha(x0,xf,h,a,E_match,m, V) 
        plt.plot(x_well_array, psi_array_well, '.', ms = 0.5)
        plt.plot(x_out_array, psi_array_out, '.', ms = 0.5)
        plt.show()
        
        print('matched, E = {0}'.format(E_match))
        
        #print(E1)
        break