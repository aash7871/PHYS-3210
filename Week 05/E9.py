#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 09:58:26 2019

@author: amandaash
"""
import numpy as np 
import matplotlib.pyplot as plt
import scipy.integrate as integral

def integrate(function, x0, xf, dx, right = True, left = False):
    x = np.arange(x0,xf,dx)
    y = function(x)
    if right == True:
        for height in y[1:]:
            summation = np.sum(y[1:])
            return dx*summation
    if left == True:
        for height in y[:-1]:
            summation = np.sum(y[:-1])
            return dx*summation
def x2(x):
    y = x**2
    return y
print('numeric integrats of $x^2$')
numeric_integral = integrate(x2, 0, 10, 0.001, True, False)
print('analytic solution = {0}'.format((10**3)/3))
print('numerical solution = {0}'.format(numeric_integral))
trap = np.trapz((np.arange(0,10,0.001))**2, np.arange(0,10,0.001))
print('numpy trapezoidal solution = {0}'.format(trap))

simpson = integral.simps((np.arange(0,10,0.001))**2, np.arange(0,10,0.001))
print('scipy simpsons rule = {0}'.format(simpson))
scipy_trap = integral.trapz((np.arange(0,10,0.001))**2, np.arange(0,10,0.001))
print('scipy trapezoidal rile = {0}'.format(scipy_trap))
#romberg = integral.romb((np.arange(0,10,0.001))**2, np.arange(0,10,0.001))
#print('scipy romberg = {0}'.format(romberg))
gaussian = integral.fixed_quad(x2, 0, 10)
print('scipy fixed quadrature = {0}'.format(gaussian[0]))


def f1(x):
    y = np.sin(100*x)
    return y
print('numeric integration of sin(100*x)')
numeric_integral = integrate(f1, 0, 2*np.pi, 0.0001, True, False)
print('analytic solution = {0}'.format(0))
print('numerical solution = {0}'.format(numeric_integral))
trap = np.trapz(f1(np.arange(0,2*np.pi,0.0001)), np.arange(0,2*np.pi,0.0001))
print('numpy trapezoidal solution = {0}'.format(trap))

simpson = integral.simps(f1(np.arange(0,2*np.pi,0.0001)), np.arange(0,2*np.pi,0.0001))
print('scipy simpsons rule = {0}'.format(simpson))
scipy_trap = integral.trapz(f1(np.arange(0,2*np.pi,0.0001)), np.arange(0,2*np.pi,0.0001))
print('scipy trapezoidal rile = {0}'.format(scipy_trap))
#romberg = integral.romb((np.arange(0,10,0.001))**2, np.arange(0,10,0.001))
#print('scipy romberg = {0}'.format(romberg))
gaussian = integral.fixed_quad(f1, 0, 2*np.pi)
print('scipy fixed quadrature = {0}'.format(gaussian[0]))

plt.plot(np.arange(0,2*np.pi,0.0001), f1(np.arange(0,2*np.pi,0.0001)), '-')
plt.show()