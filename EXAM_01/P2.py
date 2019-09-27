#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 21:08:42 2019

@author: amandaash
"""
import numpy as np

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

high_low_sum = integrate(x2, 10**4, 0, -0.00001)
low_high_sum = integrate(x2, 0, 10**4, 0.00001)

print(high_low_sum, low_high_sum)
