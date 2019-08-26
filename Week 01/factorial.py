#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 11:03:03 2019

@author: amandaash
"""

def factorial(n):
    factor = 1
    for n in range(1,((2*n)+1)+1):
        print(n)
        factorial = factor*(n)
        factor = factorial
    return factorial