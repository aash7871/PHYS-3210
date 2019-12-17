#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 09:59:40 2019

@author: amandaash
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import scipy.optimize as opt

E_i = np.arange(0,225,25)
g_Ei = np.array([10.6,16.0,45.0,83.5,52.8,19.9,10.8,8.25,4.7])
sigma_i = np.array([9.34,17.9,41.5,85.5,51.5,21.5,10.8,6.29,4.14])

print(np.mean(sigma_i))

lagrange = interp.lagrange(E_i, g_Ei)
lagrange_g_Ei = lagrange(np.arange(0,205,5))

linear = interp.interp1d(E_i, g_Ei)
lin_g_Ei = linear(np.arange(0,205,5))

cubic_spline = interp.CubicSpline(E_i, g_Ei)
cubic_g_Ei = cubic_spline(np.arange(0,205,5))

plt.title('Total lagrange polynomial')
plt.plot(E_i, g_Ei, '.')
plt.plot(np.arange(0,205,5), lagrange_g_Ei, '-')
plt.savefig('lagrange.pdf')
plt.show()

print('peak = {0}'.format(np.max(lagrange_g_Ei)))

coeff = lagrange.c

def half_max(x):
    y = coeff[0]*x**8 + coeff[1]*x**7 + coeff[2]*x**6 + coeff[3]*x**5 + coeff[4]*x**4 + coeff[5]*x**3 + coeff[6]*x**2 + coeff[7]*x + coeff[8] - (np.max(lagrange_g_Ei)/2)
    
    return y

root1 = opt.newton(half_max, 50)

root2 = opt.newton(half_max, 110)

print("FWHM = {0}".format(root2 - root1))

plt.title('linear interpolation')
plt.plot(E_i, g_Ei, '.')
plt.plot(np.arange(0,205,5), lin_g_Ei, '-')
plt.savefig('linear.pdf')
plt.show()

plt.title('cubic spline interpolation')
plt.plot(E_i, g_Ei, '.')
plt.plot(np.arange(0,205,5), cubic_g_Ei, '-')
plt.savefig('cubic_spline.pdf')
plt.show()

E_i1 = E_i[0:3]
E_i2 = E_i[2:6]
E_i3 = E_i[5:]
g_Ei1 = g_Ei[0:3]
g_Ei2 = g_Ei[2:6]
g_Ei3 = g_Ei[5:]

E_i = [[E_i1], [E_i2], [E_i3]]
g_Ei = [[g_Ei1], [g_Ei2], [g_Ei3]]

for n in range(3):
    
    interp_funcn = interp.lagrange(E_i[n][0], g_Ei[n][0])
    plt.plot(np.arange(np.min(E_i[n][0]),np.max(E_i[n][0])+5, 5), interp_funcn(np.arange(np.min(E_i[n][0]),np.max(E_i[n][0])+5, 5)), '-')
    plt.plot(E_i[n][0], g_Ei[n][0], '.')
plt.title('bisected lagrange polynomial')
plt.savefig('bisected_lagrange.pdf')
plt.show()