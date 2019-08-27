# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# -*- coding: utf-8 -*-
"""
Exercise 03: Chapter 03, Kinder & Nelson

A common way to determine the value of a function is to sum over a series. 
For example, the Maclaurin series for sin(x) is

    sin(x) = x - x**3/3! + x**5/5! - x**7/7! + ...

Perform a series expansion to derive the equation above. Next, write down 
a general expression for the sum of the series that is valid between n = 0 
and n = N, where N ≥ 0. This will serve as your algorithm for summing 
the series.

One problem with the algorithm is that we do not know which value 
of N is suitable when calcualting the series. Instead of guessing, have 
your code proceed with the summation until the Nth term contributes a 
negligible amount to the final summation, say 1 part in 10**8. We call 
this _numerical convergence_.

Before writing any lines of code, discuss an approach with your neighbor 
and write out on paper how your code should proceed. Code up your approach 
in Spyder once you're done. 

Here are your tasks:

   1. Perform a Maclaurin series expansion of the function sin(_x_) to 
      derive the equation above. 
   2. Derive a generalized, finite summation form for the series based 
      on your Maclaurin series expansion.
   3. Discuss with your neighbor about how to approach coding the problem
      and write out on paper how you code should proceed. 
   4. Code your approach in Spyder once you are finished.
   5. Show that, for small values of _x_, the algorithm converges and that
      it converges to the correct value by comparing your results to the
      value determined using NumPy's sine function.
   6. Which value for _N_ was required to reach the desired precision
      to obtain numerical convergence for small values of _x_?
   7. Steadily increase _x_ and write down the relative error between your
      calculated value for sin(_x_) and the NumPy function's value. 
   8. What do you notice about the relative error?
   9. Will there be a time when the series does not numerically converge? 
      Make a figure or two from the data you generate to support your 
      conclusion.
  10. _Challenge_ How can you modify your algorithm to be valid for any
      value of _x_?
      
   1. Perform a Maclaurin series expansion of the function sin(_x_) to 
      derive the equation above. 
   2. Derive a generalized, finite summation form for the series based 
      on your Maclaurin series expansion.
   3. Discuss with your neighbor about how to approach coding the problem
      and write out on paper how you code should proceed. 
   4. Code your approach in Spyder once you are finished.
   5. Show that, for small values of _x_, the algorithm converges and that
      it converges to the correct value by comparing your results to the
      value determined using NumPy's sine function.
   6. Which value for _N_ was required to reach the desired precision
      to obtain numerical convergence for small values of _x_?
   7. Steadily increase _x_ and write down the relative error between your
      calculated value for sin(_x_) and the NumPy function's value. 
   8. What do you notice about the relative error?
   9. Will there be a time when the series does not numerically converge? 
      Make a figure or two from the data you generate to support your 
      conclusion.
  10. _Challenge_ How can you modify your algorithm to be valid for any
      value of _x_? (Phase shift)

Created on Tue Aug 20 11:02:00 2019

@author: gafeiden
"""

import numpy as np
import matplotlib.pyplot as plt

#in a for loop:

def factorial(n):
    factor = 1
    for n in range(1,((2*n)+1)+1):
        factorial = factor*(n)
        factor = factorial
    return factorial

x_val = np.arange(-10,10,0.1)



#x_val = [0.2]
summation = 0
relative_errors = []

approx_y = []
actual_y = []

for x_v in x_val:
    print('x = ',x_v)
    print("sin(x) = ", np.sin(x_v))
    actual_y.append(np.sin(x_v))
    initial_fact = 1.0
    N = 1000
    m = x_v//(2*np.pi)
    if np.sign(x_v) == -1:
        x = x_v + (2*m*np.pi)
    if np.sign(x_v) == 1:
        x = x_v - (2*m*np.pi)
    
    summation = 0.0
    
    fact_array = np.arange(0,N+1,1)
 
    for n in fact_array:
        
        denomenator = factorial(n)
            
        iteration = (((-1.0)**(n))*(x**((2*n)+1)))/denomenator
        
        summation = summation + iteration
        
        if np.abs(iteration/summation) <= (1e-16):
            print("iteration/summation = ", np.abs(iteration/summation))
            print('sin(x) approximation complete ', n)
            break
        
    #print('sin(x) approx ', summation)
             
    
    abs_err = np.sin(x) - summation
    rel_err = abs_err/x
    relative_errors.append(rel_err)
    approx_y.append(summation)

plt.plot(x_val, relative_errors,  '.', lw=1)
plt.xlabel('x')
plt.ylabel('relative error')
plt.show()

plt.plot(x_val, actual_y, '.', lw = 1, label = 'f(x) = sin(x)', alpha = 0.3, color = 'r')
plt.plot(x_val, approx_y, '.', lw = 1, label = 'sin(x) approximation', alpha = 0.3, color = 'b')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()

plt.plot(actual_y, approx_y, '.', lw = 1)
plt.ylabel('Approximation')
plt.xlabel('Actual')
plt.show()

#plt.plot(x_val, actual_y, '.', lw = 1, label = 'f(x) = sin(x)')
#plt.plot(x_val, approx_y, '.', lw = 1, label = 'sin(x) approximation', alpha = 0.3)
#plt.xlabel('x')
#plt.ylabel('y')
#plt.xlim(-10,10)
#plt.ylim(-3,3)
#plt.legend()
#plt.show()       


