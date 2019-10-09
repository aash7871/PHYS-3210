#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 10:11:19 2019

@author: amandaash
"""

import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt


E_b = np.arange(0,10,0.001)
f_Eb = (((10-E_b)**(0.5))*np.tan((10-E_b)**(0.5))) - E_b**(0.5)

plt.plot(E_b, f_Eb)
plt.ylim(-10,10)
plt.grid()
plt.axvline(x = 8.592785275230653, color = 'k', linestyle = 'dashed')
#plt.yscale('log')
plt.savefig('binding_func_10_even.pdf')
plt.show()

def binding_energy(Eb):
    f = (((10-Eb)**(0.5))*np.tan((10-Eb)**(0.5))) - Eb**(0.5)
    return f

E_b_bisect = opt.bisect(binding_energy,8,9.5)

E_b_newton = opt.newton(binding_energy, 8.5)

E_b_brents = opt.brentq(binding_energy, 8.5, 9.5)

print('root even function: bisection method = {0}, precision = {1}'.format(E_b_bisect, binding_energy(E_b_bisect)))
print('root even function: Newton method = {0}, precision = {1}'.format(E_b_newton, binding_energy(E_b_newton)))
print('root even function: brents method = {0}, precision = {1}'.format(E_b_brents, binding_energy(E_b_brents)))
def binding_energy_beta(Eb):
    f = (np.sqrt(Eb)*(np.tan(np.sqrt(10-Eb))**-1))-np.sqrt(10-Eb)
    return(f)

    
energies = np.arange(0,10,0.001)
fE = binding_energy_beta(energies)

plt.plot(energies, fE)
plt.ylim(-20,10)
plt.axvline(x = 8.592785275230653, color = 'k', linestyle = 'dashed')
plt.grid()
plt.savefig('binding_func_10_odd.pdf')
plt.show()

Eb_bisect_beta = opt.bisect(binding_energy_beta,8,9.5)
print("root odd function: bisection method = {0}".format(Eb_bisect_beta))
E_b_newton_beta = opt.newton(binding_energy_beta, 8.5)
print("root odd function: Newton method = {0}".format(E_b_newton_beta))

def binding_energy_gamma(E, A):
    f = (((A-E)**(0.5))*np.tan((A-E)**(0.5))) - E**(0.5)
    return f

energy_gamma1 = np.arange(0,10,0.001)
energy_gamma2 = np.arange(0,20,0.001)
energy_gamma3 = np.arange(0,30,0.001)

plt.plot(energy_gamma1, binding_energy_gamma(energy_gamma1, 10))
plt.ylim(-10,10)
plt.grid()
plt.show()
plt.title('Binding energy = 20')
plt.plot(energy_gamma2, binding_energy_gamma(energy_gamma2, 20))
plt.ylim(-20,20)
plt.grid()
plt.savefig('binding_func_20_even.pdf')
plt.show()
def binding_energy_gamma20(E):
    f = (((20-E)**(0.5))*np.tan((20-E)**(0.5))) - E**(0.5)
    return f
Eb_bisect_gamma20_root1 = opt.bisect(binding_energy_gamma20,5,7.5)
Eb_bisect_gamma20_root2 = opt.bisect(binding_energy_gamma20,17.7,19)

print(Eb_bisect_gamma20_root1, binding_energy_gamma20(Eb_bisect_gamma20_root1))
print(Eb_bisect_gamma20_root2, binding_energy_gamma20(Eb_bisect_gamma20_root2))
plt.title('Binding energy = 30')
plt.plot(energy_gamma3, binding_energy_gamma(energy_gamma3, 30))
plt.ylim(-30,30)
plt.grid()
plt.savefig('binding_func_30_even.pdf')
plt.show()
def binding_energy_gamma30(E):
    f = (((30-E)**(0.5))*np.tan((30-E)**(0.5))) - E**(0.5)
    return f
Eb_bisect_gamma30_root1 = opt.bisect(binding_energy_gamma30,14,16)
Eb_bisect_gamma30_root2 = opt.bisect(binding_energy_gamma30,28,29.5)

print(Eb_bisect_gamma30_root1, binding_energy_gamma30(Eb_bisect_gamma30_root1))
print(Eb_bisect_gamma30_root2, binding_energy_gamma30(Eb_bisect_gamma30_root2))

