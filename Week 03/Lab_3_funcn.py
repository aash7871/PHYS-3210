#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 11:57:14 2019

@author: amandaash
"""
def rand_protien(N, P, epsilon):
    import random as rand
    import numpy as np
    import matplotlib.pyplot as plt

    x=0
    y=0
    coordinates = [[0,0]]
    polar_coordinates = []
    non_polar_coordinates = []
    non_polar_pairs = []

    f = []
    for iteration in range(N):
        x_add = rand.randint(-1,1)
        y_add = rand.randint(-1,1)
   
        if abs(x_add) == abs(y_add):
        
            continue
    
        x_n = x + x_add
        y_n = y + y_add
    
    
        coordinate_new = [x_n, y_n]
    
        if coordinate_new in coordinates:
            continue
    
    
        else:
            coordinates.append(coordinate_new)
            probability = rand.random()
        
            if probability >= P:
                polar_coordinates.append(coordinate_new)
            elif probability <= P:
                non_polar_coordinates.append(coordinate_new)
                neighbor_pts = [[x_n+1, y_n], [x_n-1, y_n], [x_n, y_n+1], [x_n, y_n-1]]
                for point in neighbor_pts:
                    if point in non_polar_coordinates and point != coordinate:
                        f.append(1)
                        non_polar_pairs.append(point)
        
        if [x_n + 1, y_n] in coordinates and [x_n - 1, y_n] in coordinates and [x_n, y_n + 1] in coordinates and [x_n, y_n -1] in coordinates:
            break
    
        x = x_n
        y = y_n
        coordinate = [x_n, y_n]
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

    terminal_x = x_plot[-1]
    terminal_y = y_plot[-1]
    plt.plot(x_plot, y_plot, '-')
    plt.plot(polar_coordinates[:,0], polar_coordinates[:,1], '.', color = 'orange', label = 'polar monomers')
    plt.plot(non_polar_coordinates[:,0], non_polar_coordinates[:,1], '.', color = 'grey', label = 'non-polar momomers')
    #plt.plot(non_polar_pairs[:,0], non_polar_pairs[:,1], '.', color = 'purple')
    plt.legend()
    #plt.show()
    
    return chain_length, len(f)*epsilon

l = []
E = []
protien_number = []
for protien in range(500):
    length, Energy = rand_protien(10000, 0.70, 1)
    plt.title(protien)
    plt.savefig('/Users/amandaash/Desktop/PHYS_3210/Week 04/Protien_plots/{0}.png'.format(protien))
    plt.close()
    protien_number.append(protien)
    E.append(Energy)
    l.append(length)
    
np.savetxt('/Users/amandaash/Desktop/PHYS_3210/Week 04/protien_energy.txt', np.column_stack((protien_number, l, E)), header = 'protien_number' 'length of protien' 'energy of protien')
plt.plot(l, E, '.')
plt.xlabel('Protien Length')
plt.ylabel('Protien Energy')
plt.savefig('/Users/amandaash/Desktop/PHYS_3210/Week 04/Protien_energy_length.pdf')