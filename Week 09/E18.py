#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 10:00:28 2019

@author: amandaash
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

exp = np.genfromtxt("/Users/amandaash/Desktop/PHYS_3210/Week 09/pi_meson_decays.dat")


plt.plot(exp[:,0], exp[:,1], '.')
plt.bar(exp[:,0], exp[:,1])

plt.show()

uncertainty = np.sqrt(exp[:,1])


data = np.column_stack((exp, uncertainty))


plt.errorbar(data[:,0], data[:,1], yerr = data[:,2], fmt = '.')
plt.ylabel("N(t)")
plt.xlabel("t")
plt.show()

log_uncertainties = data[:,2]*np.log(data[:,1])

log_data = np.column_stack((data[:,0], np.log(data[:,1]), log_uncertainties))

epsilon = 10**-10
log_data[9,1] = epsilon
log_data[9,2] = 2
log_data[11,2] = 2


plt.errorbar(log_data[:,0], log_data[:,1], yerr = log_data[:,2], fmt = '.')
plt.ylabel("ln(N(t))")
plt.xlabel("t")


def fit_funcn(t, N0, C):
    lnNt = np.log(N0) - (C*t)
    return lnNt
print(log_uncertainties)
parameters, parameters_covariance =  opt.curve_fit(fit_funcn, log_data[:,0], log_data[:,1], sigma = log_data[:,2])
print(parameters_covariance)
#print(N0_fit, C_fit)
plt.text(20,20,"$\\tau = {0}$ ns".format(1/C))
plt.plot(log_data[:,0], fit_funcn(log_data[:,0], parameters[0], parameters[1]))
plt.savefig("linear_fit.pdf")
plt.show()

C = parameters[1]
print('Lifetime = {0} ns'.format(1/C))

def exp_funcn(t, N0, tau):
    N_t = (-N0/tau)*np.exp(-t/tau)
    return N_t
data[9,2] = 2
parameters_alpha, parameters_vari_alhpa =  opt.curve_fit(exp_funcn, data[:,0], data[:,1], sigma = data[:,2])
print(parameters_vari_alhpa)
plt.errorbar(data[:,0], data[:,1], yerr = data[:,2], fmt = '.')
plt.plot(data[:,0], exp_funcn(data[:,0], parameters_alpha[0], parameters_alpha[1]))
plt.xlabel('t')
plt.ylabel('N(t)')
plt.text(60,25, "$\\tau = {0}$ ns".format(parameters_alpha[1]))
plt.savefig("exp_fit.pdf")
plt.show()

print("Lifetime = {0}ns".format(parameters_alpha[1]))
