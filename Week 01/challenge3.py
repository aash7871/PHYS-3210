#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 20:14:35 2019

@author: amandaash
"""
import numpy as np
import matplotlib.pyplot as plt


x_val = np.arange(-10,10,0.1)



#x_val = [0.2]
summation = 0
relative_errors = []

approx_y = []
actual_y = []

for x in x_val:
    print('x = ',x)
    print("sin(x) = ", np.sin(x))
    actual_y.append(np.sin(x))
    initial_fact = 1.0
    N = 5
    #if x == 0.0:
        #continue
    
    summation = 0.0
    
    fact_array = np.arange(0,N+1,1)
 
    for n in fact_array:
        
        value = ((2*n) + 1)
        fact = (value)*initial_fact
        initial_fact = fact

        denomenator = fact
        
        iteration = (((-1)**(n))*((x)**(value)))/denomenator
    
        summation = summation + iteration
        
        
        
        if np.abs(iteration/summation) <= (1e-8):
            print("iteration/summation = ", np.abs(iteration/summation))
            print('sin(x) approximation complete ', n)
            break
        
    print('sin(x) approx ', summation)
             
    
    abs_err = np.sin(x) - summation
    rel_err = abs_err/x
    relative_errors.append(rel_err)
    approx_y.append(summation)

plt.plot(x_val, relative_errors,  '.', lw=1)
plt.xlabel('x')
plt.ylabel('relative error')
plt.show()

plt.plot(x_val, actual_y, '.', lw = 1, label = 'f(x) = sin(x)')
plt.plot(x_val, approx_y, '.', lw = 1, label = 'sin(x) approximation')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()

plt.plot(x_val, actual_y, '.', lw = 1, label = 'f(x) = sin(x)')
plt.plot(x_val, approx_y, '.', lw = 1, label = 'sin(x) approximation')
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(-10,10)
plt.ylim(-3,3)
plt.legend()
plt.show()     

    