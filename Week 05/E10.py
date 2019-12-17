#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 10:08:09 2019

@author: amandaash
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand

def pond_sampler(N): 
    
    circle_pts = []
    square_pts = []

    x = (rand.random(N) - 0.5)*2.0
    y = (rand.random(N) - 0.5)*2.0

    for sample in range(N):
        r = np.sqrt(x[sample]**2 + y[sample]**2)
        if r<= 1:
            circle_pts.append([x[sample], y[sample]])
        else: 
            square_pts.append([x[sample], y[sample]])
    
    N_pond = len(circle_pts)
    circle_pts = np.array(circle_pts)
    square_pts = np.array(square_pts)

    #plt.plot(circle_pts[:,0], circle_pts[:,1], '.')
    #plt.plot(square_pts[:,0], square_pts[:,1], '.', color = 'g')
    #plt.savefig('/Users/amandaash/Desktop/PHYS_3210/Week 05/Pretty_pond.pdf')
    #plt.show()

    return (N_pond/N)*4

#print('Area of Pond = {0}'.format(A_pond))
#print('Number of iterations = {0}'.format(N))

Area = []
for N in range(1000):
    Area.append(pond_sampler(1000))
plt.hist(Area, bins = 30)
plt.savefig('/Users/amandaash/Desktop/PHYS_3210/Week 05/pond_histogram.pdf')
plt.show()

print(np.median(Area))
print(np.std(Area))
print(np.std(Area)/np.sqrt(1000*1000))