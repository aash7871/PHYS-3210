# -*- coding: utf-8 -*-
"""
Exercise 02: Chapter 01, Kinder & Nelson

    Download this script or enter each line in an IPython console. Your
    task is to complete the script. The script reads in data from a file
    (more on that in Chapter 4) and stores it in a variable called data.
    Investigate the properties of `data` to learn its variable type and
    other pertinent information, such as the amount of data.
    
    The data is from a Monte Carlo simulation used to determine the mass
    and age of a star that hosts an exoplanet (DS Tuc). The simulation 
    works by using random walkers to explore (nominally) all possible 
    combinations of mass, age, and other parameters to identify the most
    probable combination of parameters, and thus those that describe the
    true properties of the star. In this case, there are 6 different 
    parameters that are explored by the walker: mass, age, distance,
    brightness, surface temperature, and radius.
    
    The "location" of each walker at every step in its journey through 
    6-dimensional parameter space (chain) is preserved in the flattened 
    set of values that we loaded, effectively linking each chain together,
    end to end. However, this isn't very useful. Most of the time, we 
    want to see how each walker moved on its journey and whether all of 
    the walkers converged on the same location. To do this, we need to 
    reconstruct the individual chains from the flattened chain.
   
    So, your tasks are to:
    
       1. reconstruct the set of 400 walkers that each took 5,000 steps
          through 6D parameter space and save it to an array.
          
       2. print out the shape of your array.
       
       3. evaluate whether or not the set of chains converged. Be sure
          to confirm with your professor that you got the correct plot.
          
       4. estimate the star's mass. Compare with the mass quoted in the
          paper. (https://arxiv.org/abs/1906.10703)


Created on Tue Aug 06 14:33:00 2019

@author: gafeiden
"""
import numpy as np
import matplotlib.pyplot as plt

# read in flattened data from a file and store in variable `data`.
data = np.loadtxt("/Users/amandaash/Desktop/PHYS_3210/week 01/data/DS_Tuc_A_W0400_N5000_B0000_sDart.dat.gz")

mass, age, distance, brightness, surface_temp, radius = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]
# YOUR TASK: unflatten the array to recover the original set of 400 
#            walkers each with 5,000 steps (entries). store the results
#            in a new array called `chains`.

parameters = np.array([[mass], [age], [distance], [brightness], [surface_temp], [radius]])

mass_walk = np.reshape(mass, (400,5000))
chains = mass_walk

age_walk = np.reshape(age, (400,5000))
distance_walk = np.reshape(distance, (400,5000))
brightness_walk = np.reshape(brightness, (400,5000))
temp_walk = np.reshape(surface_temp, (400,5000))
radius_walk = np.reshape(radius, (400,5000))



# YOUR TASK: print the shape of the array to the screen
print(chains.shape)

# plot each individual chain (walker that takes N steps)
for chain in chains:
    plt.plot(chain[:], '-', lw=1)
plt.show()

# YOUR TASK: find the star's predicted mass by finding median of the 
#            first column of values (index 0) from the flattened data.

star_mass = np.median(mass)
print('The star\'s mass from the median is: ' ,star_mass)

#median_chains = []
#for chain in chains:
#    median_chains.append(np.median(chain))
#star_mass_all = np.median(median_chains)
#print('The star\'s mass from all chains is: ', star_mass_all)


# CHALLENGE: find the best mass by flattening *only* the final step from
#            each chain.

best_mass = chains[:,-1].flatten()
print('The best mass from each chain is:', best_mass)

