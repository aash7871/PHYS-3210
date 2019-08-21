# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

"""
Exercise: Chapter 01, Kinder & Nelson
    Download and run this script, or enter each line into an IPython console. 
    There are a few errors in the script. Try to find and correct them. Make 
    note of each error and explain (1) why it's an error, and (2) how to fix 
    it. 
    
    NOTE: there may be more than one way to fix an error
    
    Once you have the script working, explain why the value for c has a decimal
    and the values for a and b do not.
Created on Mon Jul 22 11:02:19 2019
@author: gafeiden
"""
import numpy as np
import matplotlib.pyplot as plt

# compute the hypotenuse of a right triangle
a = 3
b = 4 # variable b has not been defined yet, == prints a boolean 
c = np.sqrt(a**2 + b**2) # numpy has the sqrt command 

# print out the result
print("The triangle's hypotenuse is:", c ) # need parenthesis around print statement

# a and b are integer types, where as np.sqrt outputs a float. 


# All bugs : (1)b == 4 is a true false statement, and b had not been defined yet, (2) the square root 
#command comes from the numpy package, (3) the print statement needs parenthesis