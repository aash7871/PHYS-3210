#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 10:08:42 2019

@author: amandaash
"""

import numpy as np
import numpy.linalg as lin
#given
W1 = 10
W2 = 20
L1 = 3
L2 = 4
L3 = 4
L = 8

#initial guess
sintheta1 = 1
sintheta2 = 1
sintheta3 = 1
costheta1 = 1
costheta2 = 1
costheta3 = 1
T1 = 1
T2 = 1
T3 = 1

f1 = (L1*costheta1)+(L2*costheta2)+(L3*costheta3)-L
f2 = (L1*sintheta1)+(L2*sintheta2)+(-L3*sintheta3)
f3 = (sintheta1**2)+(costheta1**2) -1
f4 = (sintheta2**2)+(costheta2**2) -1 
f5 = (sintheta3**2)+(costheta3**2) -1 
f6 = (-T1*costheta1) + (T2*costheta2) 
f7 = (-T1*sintheta1) + (T2*sintheta2) + W1
f8 = (-T2*costheta2)+(T3*costheta3)
f9 = (-T2*sintheta2)-(T3*sintheta3)+W2

f_vector = np.array([[f1],[f2],[f3],[f4],[f5],[f6],[f7],[f8],[f9]])

jacobian = np.array([[0,0,0,0,L1,0,L2,0,L3],[0,0,0,L1,0,L2,0,-L3,0],\
                     [0,0,0,2*sintheta1,2*costheta1,0,0,0,0],\
                     [0,0,0,0,0,2*sintheta2,2*costheta2,0,0],\
                     [0,0,0,0,0,0,0,2*sintheta3,2*costheta3],\
                     [-costheta1,costheta2,0,0,-T1,0,T2,0,0],\
                     [-sintheta1,sintheta2,0,-T1,0,T2,0,0,0],\
                     [0,-costheta2,costheta3,0,0,0,-T2,0,T3],\
                     [0,-sintheta2,-sintheta3,0,0,-T2,0,-T3,0]])


N = 1000
for iteration in range(N):
    dx = (-lin.inv(jacobian)).dot(f_vector)
    dT1 = float(dx[0])
    dT2 = float(dx[1])
    dT3 = float(dx[2])
    d_sintheta1 = float(dx[3])
    d_costheta1 = float(dx[4])
    d_sintheta2 = float(dx[5])
    d_costheta2 = float(dx[6])
    d_sintheta3 = float(dx[7])
    d_costheta3 = float(dx[8])
    
    T1 = T1 + dT1
    T2 = T2 + dT2
    T3 = T3 + dT3
    sintheta1 = sintheta1 + d_sintheta1
    costheta1 = costheta1 + d_costheta1
    sintheta2 = sintheta2 + d_sintheta2
    costheta2 = costheta2 + d_costheta2
    sintheta3 = sintheta3 + d_sintheta3
    costheta3 = costheta3 + d_costheta3

    f1 = (L1*costheta1)+(L2*costheta2)+(L3*costheta3)-L
    f2 = (L1*sintheta1)+(L2*sintheta2)+(-L3*sintheta3)
    f3 = (sintheta1**2)+(costheta1**2) -1
    f4 = (sintheta2**2)+(costheta2**2) -1 
    f5 = (sintheta3**2)+(costheta3**2) -1 
    f6 = (-T1*costheta1) + (T2*costheta2) 
    f7 = (-T1*sintheta1) + (T2*sintheta2) + W1
    f8 = (-T2*costheta2)+(T3*costheta3)
    f9 = (-T2*sintheta2)-(T3*sintheta3)+W2

    f_vector = np.array([[f1],[f2],[f3],[f4],[f5],[f6],[f7],[f8],[f9]])


    jacobian = np.array([[0,0,0,0,L1,0,L2,0,L3],[0,0,0,L1,0,L2,0,-L3,0],\
                     [0,0,0,2*sintheta1,2*costheta1,0,0,0,0],\
                     [0,0,0,0,0,2*sintheta2,2*costheta2,0,0],\
                     [0,0,0,0,0,0,0,2*sintheta3,2*costheta3],\
                     [-costheta1,costheta2,0,0,-T1,0,T2,0,0],\
                     [-sintheta1,sintheta2,0,-T1,0,T2,0,0,0],\
                     [0,-costheta2,costheta3,0,0,0,-T2,0,T3],\
                     [0,-sintheta2,-sintheta3,0,0,-T2,0,-T3,0]])
  

    tolerance = 10**-6
    
    if dT1 <= tolerance and dT2 <= tolerance\
        and d_sintheta1 <= tolerance and d_costheta1 <= tolerance\
        and d_sintheta2 <= tolerance and d_costheta2 <= tolerance\
        and d_sintheta3 <= tolerance and d_costheta3<= tolerance:
        theta1 = np.arctan(sintheta1/costheta1)*(180/np.pi)
        theta2 = np.arctan(sintheta2/costheta2)*(180/np.pi)
        theta3 = np.arctan(sintheta3/costheta3)*(180/np.pi)
        print(iteration)
        print("T1 = {0}, T2 = {1}, T3 = {2}, theta1 = {3}, theta2 = {4}, theta3 = {5}".format(T1, T2, T3, theta1, theta2, theta3))
        break


    
    