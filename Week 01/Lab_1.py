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


#A = 160000
#B = 10
#alpha = 0.45
#beta = 0.01

#vl = A * np.exp(-alpha*time) + B * np.exp(-beta*time)

#plt.plot(HIV_data[:,0], HIV_data[:,1], '.', label = 'experimental data')
#plt.plot(time, vl, '.', label = 'model')
#plt.suptitle("A = {0}, B = {1}, alpha = {2}, beta = {3}".format(A,B,alpha,beta))
#plt.legend()
#plt.show()

#pt. c from first computer lab - Physical Modelling: The latency period of HIV is ~ ten years. 1/alpha, is equal to 
#2.22, or ~1/5 the latency period of HIV. 

#1. What were the key assumptions made in the derivation and in your approach to solving the problem that 
#made the problem tractable (easy to control)? 

#A: We assume that an anti-viral medication shuts down further infection, causing only the infected T cells to matter
# in the derivation. We also assume that the rate that T cells are cleared at is given by some probability(k)*time(t).
#This gives us an equation for the rate at which T cells are cleared from the system, a greater population of T cells
#gives a faster clearance rate. Further this equation is a differential equation, which we can solve, and yields an
#exponential. We also find the rate at which virons are cleared from the system. This equation has to consider the 
#that the rate at which virons are cleared depends on both the amount of virons and the population of infected T 
#cells. Another assumption which is made is that before the anti-viral is administered, before t = 0, the rate at 
#which virons is produced is ~ equal to the rate at which they are cleared.  

#2. What would be the consequences of relaxing the assumptions listed in question 1? How might your approach 
#to solving the problem change?

#A: If the rate at which Infected T cells were cleared was not proportional to the amount present at a given time, the 
#differential equation would be different, yielding a different equation for the amount of infected T cells at a 
#given time. This in turn would affect the viron clearance rate and the number of virons present at a given time. 


#3. How did the two limiting cases help simplify the problem?

#A: The limiting cases demonstrate the behavior of the system when the rate at which infected T cells are cleared is
#very large (K_I>>K_V), and the case where the rate at which virons are cleared is very large (K_V>>K_I). In the 
#first case, this means that new virons are not being created, and the number of virons present at a given time is 
#an exponential decay function. In the second case, it turns out that the rate at which virons are cleared and 
#produced are both very large. This gives an equation for the number of virons present which is proportional 
#to the rate at which T cells are cleared. Both of these equations provide some solution to the number of virons 
#present at a given time. We can guess a trial solution where the total number of virons present at a time is the
#summation of these two formulas. Checking this with the initial differential equation, we find that this is the 
#solution to the equation. 