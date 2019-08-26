#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 25 11:13:58 2019

@author: amandaash
"""

import numpy as np
import matplotlib.pyplot as plt
#model parameters:
clearance_rate_Tcell = 0.5
initial_viral_load = 160000
clearance_rate_virons = -0.5
X = 10
time = np.arange(0,7.2,0.1)
#defining our model
viral_load = X*np.exp(-(clearance_rate_Tcell)*time) + (initial_viral_load - X)*np.exp(clearance_rate_virons*time)

HIV_data = np.load('/Users/amandaash/Desktop/PHYS_3210/Week 01/data/HIVseries.npy')

plt.plot(HIV_data[:,0], HIV_data[:,1], '.', label = 'experimental data')
plt.plot(time, viral_load, '.', alpha = 0.5, label = 'model')
plt.suptitle("initial viral load = {0}, clearance rate T cells = {1},\nclearance rate virons = {2}, X = {3}".format(initial_viral_load, clearance_rate_Tcell, clearance_rate_virons, X))
plt.legend()
plt.show()


A = 160000
B = 10
alpha = 0.5
beta = 0.45

vl = A * np.exp(-alpha*time) + B * np.exp(-beta*time)

plt.plot(HIV_data[:,0], HIV_data[:,1], '.', label = 'experimental data')
plt.plot(time, vl, '.', label = 'model')
plt.suptitle("A = {0}, B = {1}, alpha = {2}, beta = {3}".format(A,B,alpha,beta))
plt.legend()
plt.show()