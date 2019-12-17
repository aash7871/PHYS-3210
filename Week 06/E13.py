#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 09:25:38 2019

@author: amandaash
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand
import math

def x2(x):
    y = x**2
    return y

def sinx(x):
    y = np.sin(x)
    return y

def g(x):
    y = 10*(np.cos(0.25*x))**3 
    return y

def probability_dist_integral(function, N, a, b):
    accepted= []
    rejected = []
    
    #A_square = (np.abs(a) + np.abs(b))*(np.abs(function(a)) + np.abs(function(b)))
    A_square = (np.abs(b-a))*(np.abs(function(b)-function(a)))
    print(A_square)
    x_sample = rand.uniform(a, b, N)
    y_sample = rand.uniform(np.min(function(x_sample)), np.max(function(x_sample)), N)
    points = np.column_stack((x_sample, y_sample))
    
    for coordinate in points:
        if coordinate[1] > 0:
            if coordinate[1] > 0 and coordinate[1] <= function(coordinate[0]):
                accepted.append(coordinate)
            else:
                rejected.append(coordinate)
        else:
            if coordinate[1] < 0 and coordinate[1] >= function(coordinate[0]):
                accepted.append(coordinate)
            else:
                rejected.append(coordinate)

    N_accepted = len(accepted)  
    
    accepted = np.array(accepted)
    rejected = np.array(rejected)
    plt.plot(accepted[:,0], accepted[:,1], '.')
    plt.plot(rejected[:,0], rejected[:,1], '.')
    plt.show()
    
    numeric_area =  (N_accepted/N)*A_square
    return accepted, rejected, numeric_area

plt.title('$x^2$')
accepted1, rejected1, area1 = probability_dist_integral(x2, 100000, 0, 10)
print('I(x^2)~{0}'.format(area1))
plt.title('sin(x)')
accepted2, rejected2, area2 = probability_dist_integral(sinx, 100000, 0, 2*np.pi)
print('I(sin(x))~{0}'.format(area2))
plt.title('$(10*cos^2(\\frac{x}{4}))^3$')
accepted3, rejected3, area3 = probability_dist_integral(g, 100000, 1, 50)
print('I((10*cos^2(0.25*x))^3)~{0})'.format(area3))


