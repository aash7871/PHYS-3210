#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 10:28:57 2019

@author: amandaash
"""
import numpy as np
import numpy.linalg as lin

matrix_A = np.array([[4,-2,1],[3,6,-4],[2,1,8]])
inverse_A = lin.inv(matrix_A)
print(np.array([[52,17,2],[-32,30,19],[-9,8,30]])/263)
print('---------')
print(inverse_A)
print('---------')
identity = np.dot(matrix_A, inverse_A)
print(identity)
print('-------')

b1 = np.array([12,-25,32])
solution_1 = lin.tensorsolve(matrix_A, b1)
b2 = np.array([4,-10,22])
solution_2 = lin.tensorsolve(matrix_A, b2)
b3 = np.array([20,-30,32])
solution_3 = lin.tensorsolve(matrix_A, b3)

print(solution_1)
print('-------')
print(solution_2)
print('--------')
print(solution_3)
print('-------')

#eigen value problem 
alpha = 2
beta = -1
eigen_array = [[alpha, beta],[-beta, alpha]]
solution = lin.eig(eigen_array)
print(solution)
print('--------')
eigen_array2 = np.array([[-2,2,-3], [2,1,-6], [-1,-2,0]])
solution_2 = lin.eigvals(eigen_array2)
print(solution_2)