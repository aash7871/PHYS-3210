#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 10:12:43 2019

@author: amandaash
"""

import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
"""
t1 = np.arange(0,1,0.001)
t2 = np.arange(1, 1000, 0.1)

def reduced_B_funcn(t):
    m = np.tanh(m/t)
    return(m)

def reduced_B(t):
    f = np.tanh(m/t) - m
    return f

roots_below = []
for temp in t1:
    roots1 = opt.newton(reduced_B(temp), 8.5)
    roots_below.append(roots1)
"""

t = np.arange(0,0.999,0.001)
roots_t1 = []
for temp in t:
    def reduced_B(m):
        f = np.tanh(m/temp) - m
        return f 
    roots_t1.append(opt.newton(reduced_B, 1))
    print(temp, opt.newton(reduced_B, 1))
def reduced_B_05(m):
    f = np.tanh(m/0.5) - m
    return f 
print(opt.bisect(reduced_B_05, 0.8,1.3))
   #plt.vlines(0, 10, x = roots_t1)    
plt.plot(t, roots_t1, '.')
plt.xlabel('t')
plt.ylabel('m(t)')
plt.savefig('mvt.pdf')
plt.show()
m = np.arange(0,10,0.01)
def reduced_B(m,t):
    f = np.tanh(m/t) - m
    return f

for temp in t:
    plt.plot(m, reduced_B(m,temp), '.', label = temp)
plt.axhline(y = 0, color = 'k')
plt.xlim(0,2)
#plt.legend()
plt.show()

t = np.arange(0.1,4,0.5)
m = np.arange(0,10,0.01)
for temp in t:
    plt.plot(m, reduced_B(m,temp), '.', label = temp)
plt.axhline(y = 0, color = 'k')
#plt.ylim(-1,0)
plt.xlim(0,1)
plt.legend()
plt.show()