#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 10:06:52 2019

@author: amandaash
"""

import random as rand
import numpy as np
import matplotlib.pyplot as plt

x=0
y=0
coordinates = [[0,0]]
polar_coordinates = []
non_polar_coordinates = []
non_polar_pairs = []
N = 10000
P = 0.75
epsilon = 1
f = []
for iteration in range(N):
    x_add = rand.randint(-1,1)
    y_add = rand.randint(-1,1)
    #print(x_add, y_add)
    
    if abs(x_add) == abs(y_add):
        #print('skip')
        continue
    
    x_n = x + x_add
    y_n = y + y_add
    
    
    coordinate_new = [x_n, y_n]
    #if [x_n+1, y_n] in coordinates and [x_n-1, y_n] in coordinates and [x_n, y_n+1] in coordinates and [x_n, y_n-1] in coordinates:
     #   print(iteration)
        #break
    
    if coordinate_new in coordinates:
        continue
    
    
    else:
        coordinates.append(coordinate_new)
        probability = rand.random()
        #print(probability)
        if probability >= P:
            polar_coordinates.append(coordinate_new)
        elif probability <= P:
            non_polar_coordinates.append(coordinate_new)
            neighbor_pts = [[x_n+1, y_n], [x_n-1, y_n], [x_n, y_n+1], [x_n, y_n-1]]
            for point in neighbor_pts:
                if point in non_polar_coordinates and point != coordinate:
                    #print(point)
                    f.append(1)
                    non_polar_pairs.append(point)
        
    if [x_n + 1, y_n] in coordinates and [x_n - 1, y_n] in coordinates and [x_n, y_n + 1] in coordinates and [x_n, y_n -1] in coordinates:
        #print(iteration)
        break
    
    x = x_n
    y = y_n
    coordinate = [x_n, y_n]
    #        if iteration >= 0:
     #           if neighbor_pts[0] == coordinates[iteration-1]:
      #              continue
       #         elif neighbor_pts[0] in non_polar_coordinates:
        #            f.append(1)
         #       elif neighbor_pts[0] in non_polar_coordinates or neighbor_pts[1] in non_polar_coordinates or neighbor_pts[2] in non_polar_coordinates or neighbor_pts[3] in non_polar_coordinates:
          #          f.append(1)
                
#energy calculation
#f = []
#for point in range(len(non_polar_coordinates)):
    #interest_pt = np.array(non_polar_coordinates)[point]
    #neighbor_pt
    #surrounding_pts = [interest_pt[0]+1, interest_pt[1]], [interest_pt[0]-1, interest_pt[1]], [interest_pt[0], interest_pt[1]+1], [interest_pt[0], interest_pt[1]-1]
    
    #if surrounding_pts[0] in compare_array: 
     #   f.append(1)
    #if surrounding_pts[1] in compare_array: 
     #   f.append(1)
    #if surrounding_pts[2] in compare_array: 
     #   f.append(1)
    #if surrounding_pts[3] in compare_array:
     #   f.append(1)
#print(len(f))

print('f = {0}'.format(len(f)))
chain_length = len(coordinates) - 1
print('protien length = {0}'.format(chain_length))
print("E = {0}".format(len(f)*epsilon))
coordinates = np.array(coordinates)
polar_coordinates = np.array(polar_coordinates)
non_polar_coordinates = np.array(non_polar_coordinates)
non_polar_pairs = np.array(non_polar_pairs)
x_plot = coordinates[:,0]
y_plot = coordinates[:,1]

#print(coordinates)
#print(x_plot)
#print(y_plot)
terminal_x = x_plot[-1]
terminal_y = y_plot[-1]
plt.plot(x_plot, y_plot, '-')
plt.plot(polar_coordinates[:,0], polar_coordinates[:,1], '.', color = 'orange', label = 'polar monomers')
plt.plot(non_polar_coordinates[:,0], non_polar_coordinates[:,1], '.', color = 'grey', label = 'non-polar momomers')
plt.plot(non_polar_pairs[:,0], non_polar_pairs[:,1], '.', color = 'purple')
#plt.plot(terminal_x, terminal_y, '.', markersize = 10, color = 'r')
#plt.plot(0,0,  '.', markersize = 10, color = 'g')
#plt.plot(terminal_x + 1, terminal_y, '.', markersize = 10, color = 'y')
#plt.plot(terminal_x - 1, terminal_y, '.', markersize = 10, color = 'y')
#plt.plot(terminal_x , terminal_y+1, '.', markersize = 10, color = 'y')
#plt.plot(terminal_x , terminal_y-1, '.', markersize = 10, color = 'y')
plt.legend()
plt.show()

#Energy calculation


