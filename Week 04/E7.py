#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 09:59:16 2019

@author: amandaash
"""

import numpy as np
import matplotlib.pyplot as plt
import rand
plt.title('test sequence')
P = rand.powerResidue(10)
print(P)
plt.title('repeated sequence')
repeat_sequence = rand.powerResidue(20, seed=None, a = 3, c=13, M = 2, to_plot = True)
print(repeat_sequence)
plt.savefig('/Users/amandaash/Desktop/PHYS_3210/Week 04/E7_repeat_sequence.pdf')
plt.show()
number_sequences = 15
for iteration in range(number_sequences):
    seed_value = rand.powerResidue(number_sequences)[iteration]
    sequence_length = 10
    sequence = rand.powerResidue(sequence_length, seed = seed_value, to_plot = False)
    plt.plot(np.arange(0,sequence_length,1), sequence, '-')
plt.savefig('/Users/amandaash/Desktop/PHYS_3210/Week 04/E7_multi_sequence.pdf')
plt.show()

#if you plot the sequences in this way, they will all have the same seed value since they were all generated at the
#same time. If you wanted different random sequences, you would have to change the seed values for each sequence. 
plt.title('power residue v. numpy.random.rand')
power_values = rand.powerResidue(20)
numpy_values = np.random.rand(20)
plt.plot(np.arange(1,21,1), power_values, label = 'power residue function')
plt.plot(np.arange(1,21,1), numpy_values, label = 'numpy random function')
plt.legend(loc = 'lower left')
plt.show()

power_points = [rand.powerResidue(20), rand.powerResidue(20)] 
numpy_points = [np.random.rand(20), np.random.rand(20)]
plt.title('power v. numpy, scatter')
plt.plot(power_points[0], power_points[1], '.')
plt.plot(numpy_points[0], numpy_points[1], '.')
plt.savefig('/Users/amandaash/Desktop/PHYS_3210/Week 04/E7_numpyvpower.pdf')
plt.show()


plt.title('correlation test')
N = 20
correlation_sequence = rand.powerResidue(N)
plt.plot(np.arange(0,N,1), correlation_sequence)
plt.show()
xi = correlation_sequence[0]
summation = 0
for n in correlation_sequence[1:]:
    summation = summation + xi*n
    xi = n

near_neighbor = summation/N

print('near neighbor correlation test = {0}'.format(near_neighbor))