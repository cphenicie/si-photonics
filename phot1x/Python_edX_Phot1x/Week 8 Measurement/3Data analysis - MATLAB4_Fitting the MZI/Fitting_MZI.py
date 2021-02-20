# Python 2.7 script
# by Lukas Chrostowski in Matlab
# by Huixu (Chandler) Deng in Python, 2017

from __future__ import print_function # make 'print' compatible in Python 2X and 3X
from scipy.io import loadmat        # used to load MATLAB data files
import matplotlib.pyplot as plt     # used to plot data
import os                           # used to process filenames
import urllib                       # used to download files from Dropbox.com
import numpy as np
import scipy.optimize as opt

# Enter the order of the polynomial for the curve fitting:
POLY_ORDER = 9

# Download the file from Dropbox.  Dropbox requires that you have a ?dl=1 at the end of the file
# Store the file in the local directory
# Enter the Dropbox URL here.  Make sure it has a =1 at the end:
url = 'https://www.dropbox.com/s/1rvjfef4jqybc12/ZiheGao_MZI2_271_Scan1.mat?dl=1'
dL = 111.915  # [micron] Path length difference in the MZI
FileName = os.path.split(os.path.splitext(url)[0]+'.mat')[1]

print (FileName)

if not(os.path.exists(FileName)):
    print('Downloading from Dropbox')
    if sys.version_info[0] < 3:
        import urllib                       # used to download files from Dropbox.com
        urllib.urlretrieve (url, FileName) # used in Python 2.x
    else:
        from urllib.request import urlretrieve  # used to download files from Dropbox.com
        urlretrieve (url, FileName) # used in Python 3.x

PORT=1 # Which Fibre array port is the output connected to?
# Load the file, from the MATLAB format.
matData = loadmat(FileName, squeeze_me=True, struct_as_record=False)

# Read the experimental data from the MATLAB file
lam = matData['scanResults'][PORT-1].Data[:,0]/1e9
amplitude = matData['scanResults'][PORT-1].Data[:,1]

# Plot the raw data:
plt.figure()
plt.plot(lam*1e6, amplitude)
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('Transmission [dB]')
plt.autoscale(enable=True, axis='x', tight=True)
plt.autoscale(enable=True, axis='y', tight=True)
plt.title('Experimental data (raw)')
plt.show()

# Curve fit data to a polynomial for baseline correction
p=np.polyfit((lam-np.mean(lam))*1e6, amplitude, POLY_ORDER)
amplitude_baseline=np.polyval(p,(lam-np.mean(lam))*1e6)

# Perform baseline correction to flatten the spectrum
# Use the curve polynomial, and subtract from original data
amplitude_corrected = amplitude - amplitude_baseline
amplitude_corrected = amplitude_corrected + np.max(amplitude_baseline) - np.max(amplitude)

plt.figure()
plt.plot (lam*1e6, amplitude_corrected)
plt.xlabel ('Wavelength [$\mu$m]')
plt.ylabel ('Transmission [dB]')
plt.autoscale(enable=True, axis='x', tight=True)
plt.autoscale(enable=True, axis='y', tight=True)
plt.title ('Experimental data (baseline corrected)')
plt.show()

# data only within the wavelength range of interest.
lam_min = np.min(lam) # Can limit the analysis to a range of wavelengths
lam_max = np.max(lam) #  if the data on the edges is noisy
lam_max = 1.57e-6
lam1=np.arange(lam_min, lam_max, np.min(np.diff(lam)))
amplitude1=np.interp(lam1, lam, amplitude_corrected)
amplitude1[np.nonzero(amplitude1==-np.inf)]=-50 # check if there are -infinity data points

plt.figure()
plt.plot(lam1*1e6, amplitude1)
plt.xlabel ('Wavelength [$\mu$m]')
plt.ylabel ('Transmission [dB]')
plt.autoscale(enable=True, axis='x', tight=True)
plt.autoscale(enable=True, axis='y', tight=True)
plt.title ('Experimental data (baseline corrected, wavelength range)')
plt.show()

#==============================================================================
# Define the MZI transfer function
# - as a Taylor expansion around the central wavelength
# - Use units of [microns] – keeps the variables closer to 1.
# - These make the curve fitting easier.
#==============================================================================
lam0 = np.mean(lam1)*1e6	
# effective index:
def neff(nx, lam):
    return 	nx[0] + nx[1]*(lam-lam0) + nx[2]*(lam-lam0)**2
# neff([2.4, -1, 0], 1.56)  % test it.
# alpha = 1e-3;  % propagation loss [micron^-1]
# complex propagation constant
def beta(nx, alpha, lam):
    return 2*np.pi*neff(nx, lam)/lam - 1j*alpha/2*np.ones(len(lam))
# beta([2.4, -1, 0], 1e-3, [1.56, 1.57]) % test it.
# MZI transfer function
def T_MZI(X, lam):
    return 10*np.log10(0.25*np.abs(1 + np.exp(-1j*beta(X[0:3], X[3], lam)*dL))**2) +X[4]
# T_MZI([2.4, -1, 0, 1e-3], [1.56, 1.57]) % test it.

#==============================================================================
# initial function for fitting
#==============================================================================
nx_init = np.array([2.4, -1, 0]) # CHANGE THE FIRST PARAMETER
alpha_init = 1e-3 # propagation loss [micron^-1]
x0 = np.hstack((nx_init, alpha_init, 0))

plt.figure()
plt.plot(lam1*1e6, amplitude1)
plt.plot(lam1*1e6, T_MZI(x0, lam1*1e6), linewidth = 3)
plt.xlabel ('Wavelength [$\mu$m]')
plt.ylabel ('Transmission [dB]')
plt.autoscale(enable=True, axis='x', tight=True)
plt.autoscale(enable=True, axis='y', tight=True)
plt.title ('MZI model (initial parameters)')
plt.show()

# Curve fit:  
def residuals(X, y, lam):
    return y - T_MZI(X, lam)

# curve fit to find expression for neff.
xfit, resnorm = opt.leastsq(residuals, x0, args=(amplitude1, lam1*1e6))
print (xfit)

r=np.corrcoef(amplitude1,T_MZI(xfit, lam1*1e6))
r2=r[0,1]**2
print (r2)

plt.figure()
plt.plot(lam1*1e6, amplitude1)
plt.plot(lam1*1e6, T_MZI(xfit, lam1*1e6), linewidth = 3)
plt.xlabel ('Wavelength [$\mu$m]')
plt.ylabel ('Transmission [dB]')
plt.autoscale(enable=True, axis='x', tight=True)
plt.autoscale(enable=True, axis='y', tight=True)
plt.title ('MZI model (fit parameters)')
plt.show()

# Check if the fit is good.  If so, find ng
if r2 >= 0.8:
    neff_fit = neff(xfit[0:3],lam1*1e6)
    dndlam = np.diff(neff_fit)/np.diff(lam1)
    dndlam = np.hstack((dndlam, dndlam[-1]))
    ng = neff_fit - lam1*dndlam
    
    # plot ng curve
    plt.figure()
    plt.plot(lam1*1e6, ng, linewidth = 4)
    plt.xlabel ('Wavelength [$\mu$m]')
    plt.ylabel ('Group index, n_g')
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.autoscale(enable=True, axis='y', tight=True)
    plt.title ('Group index (from MZI fit)')
    plt.show()
    
    # waveguide parameters at lam0
    ng0 = xfit[0] - lam0*xfit[1]
    print (ng0)
      