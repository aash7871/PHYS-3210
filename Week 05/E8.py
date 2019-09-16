#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 10:05:07 2019

@author: amandaash
"""
import numpy as np
import matplotlib.pyplot as plt 
import random as rand

#first derivatives
x = np.arange(-2*np.pi, 2*np.pi, 0.01)
f_x = np.sin(x)
f_prime_numpy = np.gradient(f_x, 0.01)
f_prime = []

for step in range(len(list(x))):
    if step >= 0:
        derivative = (f_x[step]-f_x[step-1])/(x[step]-x[step-1])
        f_prime.append(derivative)
    else:
        continue
plt.title('first derivative sin(x)')
plt.plot(x[1:], f_prime[1:], '.', label = 'forward difference', alpha = 0.2, color = 'b')
plt.plot(x, f_prime_numpy, '.', label = 'numpy', alpha = 0.2, color = 'r')
plt.plot(x, np.cos(x), '.', label = 'exact', alpha = 0.2, color = 'g')
plt.legend()
plt.savefig('first_derivative.pdf')
plt.show()

error_forward = (np.cos(x)[1:] - f_prime[1:])/np.cos(x)[1:]
plt.title('first derivative error forward difference')
plt.plot(x[1:], error_forward, '.')
plt.savefig('forward_diff_err.pdf')
plt.show()

error_numpy = (np.cos(x) - f_prime_numpy)/(np.cos(x))
plt.title('first derivative error numpy')
plt.plot(x[:-1], error_numpy[:-1], '.')
plt.savefig('grad_err.pdf')
plt.show()

#second derivatives 
f_prime_prime_numpy = np.gradient(f_prime_numpy, 0.01)
f_prime_prime = []
for step in range(len(list(x))):
    if step >= 0:
        derivative = (f_prime[step]-f_prime[step-1])/(x[step]-x[step-1])
        f_prime_prime.append(derivative)
    else:
        continue
plt.title('second derivative sin(x)')
plt.plot(x[2:], f_prime_prime[2:], '.', label = 'forward difference', alpha = 0.2, color = 'b')
plt.plot(x, f_prime_prime_numpy, '.', label = 'numpy', alpha = 0.2, color = 'r')
plt.plot(x, -np.sin(x), '.', label = 'exact', alpha = 0.2, color = 'g')
plt.legend()
plt.savefig('2nd_derivative.pdf')
plt.show()

error_forward_2 = (-np.sin(x)[2:] - f_prime_prime[2:])/(-np.sin(x)[2:])
plt.title('second derivative error forward difference')
plt.plot(x[2:], error_forward_2, '.')
plt.savefig('2nd_derivative_err_forward_diff.pdf')
plt.show()

error_numpy_2 = (-np.sin(x) - f_prime_prime_numpy)/-np.sin(x)
plt.title('second derivative error numpy')
plt.plot(x[1:-2], error_numpy_2[1:-2], '.')
plt.savefig('2nd_derivative_err_numpy.pdf')
plt.show()

#random noise
epsilon_array = []
for n in range(len(list(x))):
    epsilon_array.append(rand.random()*0.001)
f_x = np.sin(x) + epsilon_array
f_prime_numpy = np.gradient(f_x, 0.01)
f_prime = []

for step in range(len(list(x))):
    if step >= 0:
        derivative = (f_x[step]-f_x[step-1])/(x[step]-x[step-1])
        f_prime.append(derivative)
    else:
        continue
plt.title('first derivative sin(x) + $\epsilon$')
plt.plot(x[1:], f_prime[1:], '.', label = 'forward difference', alpha = 0.2, color = 'b')
plt.plot(x, f_prime_numpy, '.', label = 'numpy', alpha = 0.2, color = 'r')
plt.plot(x, np.cos(x), '.', label = 'exact', alpha = 0.2, color = 'g')
plt.legend()
plt.savefig('noisy_derivative.pdf')
plt.show()

error_forward = np.cos(x)[1:] - f_prime[1:]
plt.title('first derivative error forward difference')
plt.plot(x[1:], error_forward, '.')
plt.show()

error_numpy = np.cos(x) - f_prime_numpy
plt.title('first derivative error numpy')
plt.plot(x[:-1], error_numpy[:-1], '.')
plt.show()

#second derivatives 
f_prime_prime_numpy = np.gradient(f_prime_numpy, 0.01)
f_prime_prime = []
for step in range(len(list(x))):
    if step >= 0:
        derivative = (f_prime[step]-f_prime[step-1])/(x[step]-x[step-1])
        f_prime_prime.append(derivative)
    else:
        continue
plt.title('second derivative sin(x) + $\epsilon$')
plt.plot(x[2:], f_prime_prime[2:], '.', label = 'forward difference', alpha = 0.2, color = 'b')
plt.plot(x, f_prime_prime_numpy, '.', label = 'numpy', alpha = 0.2, color = 'r')
plt.plot(x, -np.sin(x), '.', label = 'exact', alpha = 0.2, color = 'g')
plt.legend()
plt.show()

error_forward_2 = -np.sin(x)[2:] - f_prime_prime[2:]
plt.title('second derivative error forward difference')
plt.plot(x[2:], error_forward_2, '.')
plt.show()

error_numpy_2 = -np.sin(x) - f_prime_prime_numpy
plt.title('second derivative error numpy')
plt.plot(x[1:-2], error_numpy_2[1:-2], '.')
plt.show()
