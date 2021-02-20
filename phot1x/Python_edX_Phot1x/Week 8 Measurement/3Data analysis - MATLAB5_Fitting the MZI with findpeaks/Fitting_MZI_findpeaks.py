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
from scipy.signal import savgol_filter
import sys

# Download the file from Dropbox.  Dropbox requires that you have a ?dl=1 at the end of the file
# Store the file in the local directory
# Enter the Dropbox URL here.  Make sure it has a =1 at the end:

# TE
url = 'https://www.dropbox.com/s/1rvjfef4jqybc12/ZiheGao_MZI2_271_Scan1.mat?dl=1'
dL = 111.915 # [micron] Path length difference in the MZI

# TM:
#url = 'https://www.dropbox.com/s/onjwrarapf6dumv/ADCL_1_1153_Scan1_L2-L1%3D100um.mat?dl=1'
#dL = 100

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

# Curve fit data to a polynomial for baseline correction
p=np.polyfit((lam-np.mean(lam))*1e6, amplitude, 4)
amplitude_baseline=np.polyval(p,(lam-np.mean(lam))*1e6)

# Perform baseline correction to flatten the spectrum
# Use the curve polynomial, and subtract from original data
amplitude_corrected = amplitude - amplitude_baseline
amplitude_corrected = amplitude_corrected + np.max(amplitude_baseline) - np.max(amplitude)

# data only within the wavelength range of interest.
lam_min = np.min(lam) # Can limit the analysis to a range of wavelengths
lam_max = np.max(lam) #  if the data on the edges is noisy
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
# Find peaks, extract FSR and ng, and neff
# as initial parameters
#==============================================================================
# smooth the data to find peaks accurately
amplitude_smooth=savgol_filter(amplitude1, 2001, 5)
from scipy.signal import argrelmax
indexes = argrelmax(-amplitude_smooth)

# peak amplitude above half the whole amplitude is incorrect and removed
amplitude_half = (np.max(amplitude_smooth) + np.min(amplitude_smooth))/2
index_remove = np.nonzero(amplitude_smooth[indexes] > amplitude_half)
indexes = np.delete(indexes, index_remove)
x_values = lam1[indexes]

# plot spectrum with peaks
plt.figure()
plt.plot(x_values*1e6, np.interp(x_values, lam1, amplitude1), 'ro')
plt.plot(lam1*1e6, amplitude1)
plt.xlabel ('Wavelength [$\mu$m]')
plt.ylabel ('Transmission [dB]')
plt.autoscale(enable=True, axis='x', tight=True)
plt.autoscale(enable=True, axis='y', tight=True)
plt.title ('Experimental data (with peaks)')
plt.show()

# Calculate and plot ng data points from FSR
lam_ng = (x_values[0:len(x_values)-1] + x_values[1:])/2
FSR = (x_values[1:] - x_values[0:len(x_values)-1])
ng = np.abs(lam_ng**2/1/(dL*1e-6)/FSR)

# find average ng from all reasonable ng values:
indexes = np.nonzero(ng >= 3)
ng = ng[indexes]
lam_ng = lam_ng[indexes]
FSR = FSR[indexes]

indexes = np.nonzero(ng < 5)
ng = ng[indexes]
lam_ng = lam_ng[indexes]
FSR = FSR[indexes]
ng_av = np.mean(ng)

if len(ng) == 0:
  print ('No group index values found.  Terminating script.')

print ('(estimate) Group index: ', ng_av)

# plot FSR
plt.figure()
plt.plot (lam_ng*1e6, FSR*1e9, '-o', linewidth = 1)
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('Free Spectral Range [nm]')
plt.xlim(np.array([np.min(lam1), np.max(lam1)])*1e6)
plt.title('Free Spectral Range, from Experimental data')
plt.show()

# plot ng
plt.figure()
plt.plot (lam_ng*1e6, ng, '-o',linewidth = 1)
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel ('Group Index')
plt.xlim(np.array([np.min(lam1), np.max(lam1)])*1e6)
plt.title ('Waveguide Group Index, from Experimental data')
plt.show()


#==============================================================================
# find starting point for curve fitting MZI, using the ng data
#==============================================================================
# Part 1 - n1
lam0 = x_values[int(len(x_values)/2)]*1e6 # lam0 is in microns.
n1_initial = 2.4
modeNumber = n1_initial*dL/lam0 - 0.5
n1 = (2*np.rint(modeNumber)+1)*lam0/2/dL

# Part 1 - n2 from ng_av
n2 = (n1-ng_av)/lam0

# Part 3 - n3 from slope of ng vs. lam, to get the dispersion
def f_line(x, xdata):
    return x[1]*xdata + x[0]

def residuals(x, y, xdata):
    return y - f_line(x, xdata)

[xfit, resnorm] = opt.leastsq(residuals, np.array([ng_av, 0]), args=(ng, lam_ng*1e6))
print (xfit)

r=np.corrcoef(amplitude1,f_line(xfit, lam1*1e6))
r2_ng=r[0,1]**2
print ('Goodness of fit, r^2 value: ', r2_ng)


plt.plot(lam_ng*1e6, f_line(xfit, lam_ng*1e6), linewidth = 3)
n3 = -xfit[1]/2/lam0
Dispersion0 = -(lam0*1e-6) / 299792458 * 2* (n3*1e12) * 1e12 /1e9 /1e-3 	
print ('(estimate from ng slope) Dispersion [ps/nm/km]: ', Dispersion0)

# only use the result if the fit is good
if r2_ng < 0.01:
    n3 = 0

# Initial conditions for fitting:
nx_init = np.array([n1, n2, n3])
alpha_init = 1e-3  # propagation loss [micron^-1]
x0 = np.hstack((nx_init, alpha_init, 0))

#==============================================================================
# Define the MZI transfer function
# - as a Taylor expansion around the central wavelength
# - Use units of [microns] – keeps the variables closer to 1.
# - These make the curve fitting easier.
#==============================================================================
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
    print ('Waveguide parameters at wavelength [um]: ', lam0)
    ng0 = xfit[0] - lam0*xfit[1]
    print ('Group index:', ng0)
    
    # Dispersion:
    c = 299792458
    Dispersion0 = -(lam0*1e-6)/c*2*(xfit[2]*1e12) # [s/m^2]
    # [ps/nm/km]: 
    Dispersion0 = Dispersion0 * 1e12 /1e9 /1e-3
    print ('Dispersion [ps/nm/km]: ', Dispersion0)
