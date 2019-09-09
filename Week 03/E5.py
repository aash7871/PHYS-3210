#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 09:57:04 2019

@author: amandaash
"""
import numpy as np
import matplotlib.pyplot as plt

def harmonic(N):
    """
    input: 
        N - number of iterations in the summation.
    output: 
        harmonic1 - upward harmonic series
        harmonic2 - downward harmonic series
    """
    harmonic1 = 0
    harmonic2 = 0
    for iteration in range(1,N+1,1):
        harmonic1 = harmonic1 + (1/iteration)
    for iteration in range(N,0,-1):
        harmonic2 = harmonic2 + 1/iteration
    return harmonic1, harmonic2

n_val = []
harm_err = []
for n in range(1000):
    n_val.append(n)
    harm1, harm2 = harmonic(n)
    #print(harmonic1, harmonic2)
    harmonic_error = (harm1-harm2)/(np.abs(harm1)+np.abs(harm2))
    harm_err.append(harmonic_error)

plt.plot(n_val, harm_err, '.', lw = 1)
plt.show()

"""
When N = 6, the analytical solution to the harmonic series is 2.45. The downward harmonic series is equal to 
this value, making the downward solution more precise. The upward harmonic series has extra digits. 

The reason this is the case is because when you start with a large number, in this case 1, if you are adding 
a smaller number, the smaller number is converted to scientific notation with the same exponent as the larger number.
This can truncate numbers which are greater than 16 decimal points away. For example, small numbers which would 
normally round to some number would be saved as 0. However, if you start with a small number, as in the downward 
series, the exponents are smaller initially, so more f the decimal points are saved. 

"""