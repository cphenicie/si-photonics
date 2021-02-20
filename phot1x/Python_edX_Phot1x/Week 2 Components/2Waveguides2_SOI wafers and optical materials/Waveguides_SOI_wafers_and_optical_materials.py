# Python 2.7 script
# by Lukas Chrostowski in Matlab, 2015
# by Huixu (Chandler) Deng in Python, 2017

from __future__ import print_function # make 'print' compatible in Python 2X and 3X
import matplotlib.pyplot as plt
import numpy as np
        
delta_T = 10
Thermal_coefficient = 0.01 # replace with correct value
delta_n = Thermal_coefficient * delta_T

print ('delta_n =%.8f' %delta_n)
        
      