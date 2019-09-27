#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 10:17:04 2019

@author: amandaash
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand

def MC_multidimension_integrand(function, intervals, N_dimensions, sample_number):
    """
    INPUTS:
        - function - N-dimenstional function to integrate
        - intervals - array or intervals to sample over N dimensions
        - N_dimensions - number of dimensions of function
        - number of samples taken in each dimension
    OUTPUT: 
        - numeric approximation of N-dimensional integrand
    """
    
    interval_samples = []
    interval_scalar = []
    for dimension in intervals:
        a = dimension[0]
        b = dimension [1]
        sample = rand.uniform(a, b, sample_number)
        interval_samples.append(sample)
        interval_scalar.append(b-a)
    
    interval_product = np.product(interval_scalar)
    
    interval_samples = np.array(interval_samples)
    
    interval_samples = np.column_stack((interval_samples))

    
    f_sample = []
    
    for x_sample in interval_samples:
        
        f_sample.append(function(x_sample))
        
    funcn_sum = np.sum(f_sample)
    
    numeric_integral = (interval_product/sample_number)*funcn_sum
    
    return numeric_integral

def f_ND(D_array):
    x1 = D_array[0]
    x2 = D_array[1]
    x3 = D_array[2]
    x4 = D_array[3]
    x5 = D_array[4]
    x6 = D_array[5]
    x7 = D_array[6]
    x8 = D_array[7]
    x9 = D_array[8]
    x10 = D_array[9]
    y = (x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)**2
    return y


N_interval = [[0,1], [0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1]]
D10_integral = MC_multidimension_integrand(f_ND, N_interval, 10, 10000)
print(D10_integral)

sample_array = np.arange(2,10000,10)
integral_value = []
for s in sample_array:
    num_integral = MC_multidimension_integrand(f_ND, N_interval, 10, s)
    integral_value.append(np.abs((155/6)-num_integral)/(num_integral))

plt.plot(1/np.sqrt(sample_array), integral_value, '.')
plt.xlabel(r'$\frac{1}{\sqrt{N}}$')
plt.ylabel('relative error')
plt.savefig('rel_erro.pdf')
plt.show()

    
    