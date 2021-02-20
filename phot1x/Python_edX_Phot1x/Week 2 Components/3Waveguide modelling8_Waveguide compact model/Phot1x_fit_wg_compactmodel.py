# Python 2.7 script
# by Lukas Chrostowski in Matlab
# by Huixu (Chandler) Deng in Python, 2017

#==============================================================================
# Phot1x_fit_wg_compactmodel
# by Lukas Chrostowski, 2015
# 
# User provides a matrix of neff values vs. wavelength
# curve fits to an expression.
# 
# url='https://www.dropbox.com/s/xv4he4preyfa9v2/wg-export-TM.mat?dl=1'
#==============================================================================

from __future__ import print_function # make 'print' compatible in Python 2X and 3X
from scipy.io import loadmat        # used to load MATLAB data files
import matplotlib.pyplot as plt     # used to plot data
import os                           # used to process filenames
import urllib                       # used to download files from Dropbox.com
import numpy as np
import scipy.optimize as opt

# Download the file from Dropbox.  Dropbox requires that you have a ?dl=1 at the end of the file
# Store the file in the local directory

url = "https://www.dropbox.com/s/xv4he4preyfa9v2/wg-export-TM.mat?dl=1"
FileName = os.path.split(os.path.splitext(url)[0]+'.mat')[1]
print (FileName)
urllib.urlretrieve (url, FileName) # used in Python 2.x
#urllib.request.urlretrieve (url, FileName) # used in Python 3.x

# Load the file, from the MATLAB format.
matData = loadmat(FileName, squeeze_me=True, struct_as_record=False)

neff = np.real(matData['neff']) # take the real part of the effective index
f = matData['f'] # f is the matrix of frequency points, where the effective index is recorded.
c=299792458.0  # speed of light, m/s
lams = c/f*1e6 # convert to microns
lam0 = 1.55;   # replace with desired centre wavelength

# function for the effective index expression:
def neff_eq(nx, lam):
    return nx[0] + nx[1]*(lam-lam0) + nx[2]*(lam-lam0)**2

# In Python, to do the curve fitting, 
# the leastsq function is used and the residuals between the data and the model should be created.
# function for residuals between the data and the model
def residuals(nx, y, lam):
    return y - neff_eq(nx, lam)

# initial guess
nx0=np.array([1.7,0,0])

# curve fit to find expression for neff.

nx, flag = opt.leastsq(residuals, nx0, args=(neff, lams))

r=np.corrcoef(neff,neff_eq(nx, lams))
r2=r[0,1]**2
print ('Goodness of fit, r^2 value: ', r2)
print ('n1 = ', nx[0], 'n2 =', nx[1], 'n3 =', nx[2])

lams2=np.linspace(np.min(lams), np.max(lams), 100)

plt.figure()
plt.plot(lams, neff,'o', label='Data')
plt.plot(lams, neff_eq(nx0, lams), 'r', label='Initial Guess')
plt.plot(lams2, neff_eq(nx, lams2), 'k', label='Curve Fit')
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('Effective Index')
plt.legend()
plt.show()
