#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 10:00:21 2019

@author: amandaash
"""

#step (1) Sample the function x**2 1000 and 1000 times over the interval 0 to 10
#step (2) sum the samples
#step (3) multiply sum by b-a/number of samples
#step (4) viola an integral

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand

def x2(x):
    y = x**2
    return y

def sampler(funcn, N, a, b):
    sample = rand.uniform(a, b, N)
    funcn_sample = funcn(sample)
    I = ((b-a)/N)*np.sum(funcn_sample)
    return I, np.array(sample), np.array(funcn_sample)

I_1000, x_1000, y_1000 = sampler(x2, 1000, 0, 10)
I_10000, x_10000, y_10000 = sampler(x2, 10000, 0, 10)

print('1000 samples: I ~ {0}'.format(I_1000))
print('10000 samples: I ~ {0}'.format(I_10000))

plt.plot(x_1000, y_1000, '.')
plt.title('1000 samples')
plt.savefig('1000_samples.pdf')
plt.show()
plt.plot(x_10000, y_10000, '.')
plt.title('10000 samples')
plt.savefig('100000 samples.pdf')
plt.show()

N = 10000
I_values = []
y_values = []
y_std = []
for iteration in range(N):
    I, x, y = sampler(x2, 1000, 0, 10)
    I_values.append(I)
    y_values.append(y)
    y_std.append(np.std(y))

y_values = np.array(y_values).reshape(N*1000)
print(np.std(y_values)/np.sqrt(N*1000))
print(np.std(y_std)/np.sqrt(N))
    
plt.hist(I_values, bins = 25)
plt.axvline(np.median(I_values), label = 'median I = {0}'.format(str(np.median(I_values))[:7]), color = 'r')
plt.axvline(np.mean(I_values), label = 'mean I = {0}'.format(str(np.mean(I_values))[:7]), color = 'orange')
plt.axvline(np.mean(I_values)+np.std(I_values), color = 'purple')
plt.axvline(np.mean(I_values)-np.std(I_values), color = 'purple')
plt.title('{0} samples'.format(1000*N))
plt.legend()
plt.savefig('sample hist')
plt.show()

print("{0} samples: I = {1} +/- {2}".format(N*1000, str(np.mean(I_values))[:7], str(np.std(I_values)/np.sqrt(N*1000))[:5]))


def f1(x):
    y = np.sin(100*x)
    return y
I_sin, x_sin, y_sin = sampler(f1, 10000000, 0, 2*np.pi)
print('10000 samples sin(x): I ~ {0}'.format(I_sin))
plt.plot(x_sin, y_sin, '.')
plt.show() 