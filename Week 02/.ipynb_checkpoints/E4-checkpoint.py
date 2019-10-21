#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 08:19:09 2019

@author: amandaash
"""

import numpy as np
import matplotlib.pyplot as plt

def viral_load(A,B,alpha,beta, time):
    
    V = A * np.exp(-alpha*time) + B * np.exp(-beta*time)
    
    return V

A = [160000, 3000, 50000, 20000, 100000]
B = [10, 10, 10, 10, 10]
alpha = [0.45, -0.45, -0.1, 0.1, 5]
beta = [0.01, 0.01, 0.01, 0.01, 0.01]
t = np.arange(0,7.2,0.1)

parameter_array = np.column_stack((A,B,alpha,beta))

for parameter_set in parameter_array:
    
    virus = viral_load(parameter_set[0], parameter_set[1], parameter_set[2], parameter_set[3], t)
    
    plt.plot(t, virus, '.', label = r'A = {0}, B = {1}, $\alpha$ = {2}, $\beta$ = {3}'.format(parameter_set[0], parameter_set[1], parameter_set[2], parameter_set[3]))
    

plt.ylabel('viral load')
plt.xlabel('time')
plt.title('HIV virus concentration over time')
plt.legend()
plt.tight_layout()
plt.savefig('/Users/amandaash/Desktop/PHYS_3210/Week 02/HIV_plot.pdf')
plt.show()