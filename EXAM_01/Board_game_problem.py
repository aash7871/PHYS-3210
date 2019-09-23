#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 21:21:38 2019

@author: amandaash
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand


def spin(version, N_spins):
    """
    Inputs: 
        - version - board game version, equal probability for each spin (1), half probability for deductions (2)
        - N_spins - maximum number of spins during game
    Outputs: 
        - basket - array of outcome from each spin (-all, -2, -1, 0, +1, +2, +3, +4)
        - number_cherries - array of the number of cherries in the basket at each turn
        - N - turn where game is won, sum of basket == 10
    """
    basket = []
    number_cherries = [0]
    if version == 1:
        P1 = float(1/7)
        P2 = float(2/7)
        P3 = float(3/7)
        P4 = float(4/7)
        P5 = float(5/7)
        P6 = float(6/7)
        P7 = float(1)
    if version == 2:
        P1 = float(2/11)
        P2 = float(4/11)
        P3 = float(6/11)
        P4 = float(8/11)
        P5 = float(9/11)
        P6 = float(10/11)
        P7 = float(1)
    
    for N in range(N_spins):
        P_spin = rand.random()
        if P_spin <= P1:
            basket.append(1)
        if P_spin <= P2 and P_spin > P1:
            if np.sum(np.array(basket)) == 9:
                basket.append(1)
            else:
                basket.append(2)
        if P_spin <= P3 and P_spin > P2:
            if np.sum(np.array(basket)) <= 9 and np.sum(np.array(basket)) >= 8:
                basket.append(10- np.sum(np.array(basket)))
            else:
                basket.append(3)
        if P_spin <= P4 and P_spin > P3:
            if np.sum(np.array(basket)) <= 9  and np.sum(np.array(basket)) >= 7:
                basket.append(10- np.sum(np.array(basket)))
            else:
                basket.append(4)
        if P_spin <=P5 and P_spin > P4:
            if np.sum(np.array(basket)) >= 2:
                basket.append(-2)
            if np.sum(np.array(basket)) == 1:
                basket.append(-1)
            if np.sum(np.array(basket)) == 0:
                basket.append(0)
        if P_spin <=P6 and P_spin > P5:
            if np.sum(np.array(basket)) >= 2:
                basket.append(-2)
            if np.sum(np.array(basket)) == 1:
                basket.append(-1)
            if np.sum(np.array(basket)) == 0:
                basket.append(0)
        if P_spin <= P7 and P_spin > P6:
            if np.sum(np.array(basket)) > 0:
                basket.append(-np.sum(np.array(basket)))
            else:
                basket.append(0)
        
        number_cherries.append(np.sum(np.array(basket)))
        
        if np.sum(np.array(basket)) >= 10:
            return basket, number_cherries, N
            break

plays = 10000
wins1 = []
wins2 = []
for games in range(plays):
   number1, basket1, spin_to_win1 = spin(1, 500)
   wins1.append(spin_to_win1)
   number2, basket2, spin_to_win2 = spin(2,500)
   wins2.append(spin_to_win2)

plt.hist(wins1, bins = 25, alpha = 0.3, color = '#2E2EFE', label = 'turns to win V1')
plt.hist(wins2, bins = 25, alpha = 0.3, color = '#8000FF', label = 'turns to win V2')
plt.axvline(x = np.median(wins1), color = '#2E2EFE', label = 'median spins to win V1 = {0}'.format(np.median(wins1)))
plt.axvline(x = np.median(wins2), color = '#8000FF', label = 'median spins to win V2 ={0}'.format(np.median(wins2)))
plt.xlabel('spins to win game')
plt.ylabel('number of games')
plt.legend()
plt.savefig('/Users/amandaash/Desktop/PHYS_3210/Exam_01/4b_hist.pdf')
plt.show()

plt.plot(basket1, '-', color = '#2E2EFE', label = 'example game V1')
plt.plot(basket2, '-', color = '#8000FF', label = 'example game V2')
plt.plot(basket1, '.', color = '#2E2EFE')
plt.plot(basket2, '.', color = '#8000FF')
plt.ylabel('turn outcome')
plt.xlabel('turn number')
plt.legend()
plt.savefig('/Users/amandaash/Desktop/PHYS_3210/Exam_01/4b_ex_game.pdf')
plt.show()

print('Average number of spins to win game version 1 = {0} +/- {1}, median number of spins to win V1 = {2}'.format(np.mean(wins1), str(np.std(wins1)/np.sqrt(plays))[:5], np.median(wins1)))
print('Average number of spins to win game version 2 = {0} +/- {1}, median number of spins to win V2 = {2}'.format(np.mean(wins2), str(np.std(wins2)/np.sqrt(plays))[:5], np.median(wins2)))



         
