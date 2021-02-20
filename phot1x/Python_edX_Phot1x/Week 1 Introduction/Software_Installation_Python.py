# Python 2.7 script
# by Lukas Chrostowski in Matlab, 2015
# by Huixu (Chandler) Deng in Python, 2017

from __future__ import print_function 
# make 'print' compatible in Python 2X and 3X
import matplotlib.pyplot as plt
import numpy

a = 1
b = 2
c = a + b

print ('a=', a)
print ('b=', b)
print ('c=', c)

# Practice figures:
x = numpy.arange(1,10.1,0.1)

plt.figure()
plt.plot(x, numpy.sin(x))
plt.title('The First figure')

plt.figure()
plt.plot(x, numpy.exp(x))
plt.title('The Second figure')