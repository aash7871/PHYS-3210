#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 08:19:43 2019

@author: amandaash
"""

import random as rand
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
def rand_protien_3D(N, P, epsilon, upper_lim):
    import random as rand
    import numpy as np
    import matplotlib.pyplot as plt

    x=0
    y=0
    z=0
    coordinates = [[0,0,0]]
    polar_coordinates = []
    non_polar_coordinates = []
    non_polar_pairs = []

    f = []
    for iteration in range(N):
        #x_add = rand.randint(-1,1)
        #y_add = rand.randint(-1,1)
        #z_add = rand.randint(-1,1)
        
        rand_xyz = rand.randint(0,2)
        if rand_xyz == 0:
            x_add = 1
            y_add = 0
            z_add = 0
        if rand_xyz == 1:
            x_add = 0
            y_add = 1
            z_add = 0
        if rand_xyz == 2:
            x_add = 0
            y_add = 0
            z_add = 1
   
        #if abs(x_add) == abs(y_add) or abs(x_add) == abs(z_add) or abs(y_add) == abs(z_add) or abs(x_add) == abs(z_add) == abs(y_add): 
            #continue
        
        x_n = x + x_add
        y_n = y + y_add
        z_n = z + z_add
    
        coordinate_new = [x_n, y_n, z_n]
    
        if coordinate_new in coordinates:
            continue
        #if abs(x_n) >= upper_lim or abs(y_n)>= upper_lim or abs(z_n)>= upper_lim:
         #   continue
    
        else:
            coordinates.append(coordinate_new)
            probability = rand.random()
        
            if probability >= P:
                polar_coordinates.append(coordinate_new)
            elif probability <= P:
                non_polar_coordinates.append(coordinate_new)
                neighbor_pts = [[x_n+1, y_n, z_n], [x_n-1, y_n, z_n], [x_n, y_n+1, z_n], [x_n, y_n-1, z_n], [x_n, y_n, z_n+1], [x_n, y_n, z_n-1]]
                for point in neighbor_pts:
                    if point in non_polar_coordinates and point != coordinate:
                        f.append(1)
                        non_polar_pairs.append(point)
        
        if [x_n + 1, y_n, z_n] in coordinates and [x_n - 1, y_n, z_n] in coordinates and [x_n, y_n + 1, z_n] in coordinates and [x_n, y_n -1, z_n] in coordinates and [x_n, y_n, z_n+1] in coordinates and [x_n, y_n, z_n-1] in coordinates:
            print(iteration)
            break
    
        x = x_n
        y = y_n
        z = z_n
        coordinate = [x_n, y_n, z_n]
    #print('f = {0}'.format(len(f)))
    chain_length = len(coordinates) - 1
    #print('protien length = {0}'.format(chain_length))
    #print("E = {0}".format(len(f)*epsilon))
    coordinates = np.array(coordinates)
    polar_coordinates = np.array(polar_coordinates)
    non_polar_coordinates = np.array(non_polar_coordinates)
    non_polar_pairs = np.array(non_polar_pairs)
    x_plot = coordinates[:,0]
    y_plot = coordinates[:,1]
    z_plot = coordinates[:,2]

    terminal_x = x_plot[-1]
    terminal_y = y_plot[-1]
    terminal_z = z_plot[-1]
    fig = plt.figure()
    ax = plt.subplot(111, projection='3d')
    ax.plot(x_plot, y_plot, z_plot, '-')
    ax.plot(polar_coordinates[:,0], polar_coordinates[:,1], polar_coordinates[:,2], '.', color = 'orange', label = 'polar monomers')
    ax.plot(non_polar_coordinates[:,0], non_polar_coordinates[:,1], non_polar_coordinates[:,2], '.', color = 'grey', label = 'non-polar momomers')
    #plt.plot(non_polar_pairs[:,0], non_polar_pairs[:,1], '.', color = 'purple')
    ax.legend()
    #plt.show()
    
    return chain_length, len(f)*-epsilon

l = []
E = []
protien_number = []
for protien in range(5):
    length, Energy = rand_protien_3D(1000, 0.70, 1, 10)
    plt.title(protien)
    plt.savefig('/Users/amandaash/Desktop/PHYS_3210/Week 04/protien_plots_3D/{0}.pdf'.format(protien))
    plt.close()
    protien_number.append(protien)
    E.append(Energy)
    l.append(length)
    
np.savetxt('/Users/amandaash/Desktop/PHYS_3210/Week 04/protien_energy_3D.txt', np.column_stack((protien_number, l, E)), header = 'protien_number' 'length of protien' 'energy of protien')
plt.plot(l, E, '.')
plt.xlabel('Protien Length')
plt.ylabel('Protien Energy')
plt.title('3D protein simulation')
plt.savefig('/Users/amandaash/Desktop/PHYS_3210/Week 04/Protien_energy_length_3D.pdf')