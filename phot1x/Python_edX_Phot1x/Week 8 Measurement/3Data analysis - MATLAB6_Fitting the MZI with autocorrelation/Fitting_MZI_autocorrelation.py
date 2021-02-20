# Python 2.7 script
# by Lukas Chrostowski in Matlab
# by Huixu (Chandler) Deng in Python, 2017

from __future__ import print_function # make 'print' compatible in Python 2X and 3X
from scipy.io import loadmat        # used to load MATLAB data files
import matplotlib.pyplot as plt     # used to plot data
import os                           # used to process filenames
import numpy as np
import scipy.optimize as opt
import sys
#from scipy.signal import savgol_filter

#Enter the Dropbox URL here.  Make sure it has a =1 at the end:
url = 'https://www.dropbox.com/s/1rvjfef4jqybc12/ZiheGao_MZI2_271_Scan1.mat?dl=1';
dL = 111.915 # [micron] Path length difference in the MZI


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
lam_max = 1.57e-6;
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
# Find ng from autocorrelation-based frequency estimation of spectrum
# auto-correction
#==============================================================================
r = np.correlate(amplitude1, amplitude1, 'full')
r = r[int(len(r)/2):]
lags = np.arange(len(r))

r = r[np.nonzero(lags >= 0)]
lags = lags[np.nonzero(lags >= 0)]

plt.figure()
plt.plot(lags,r)

# estimate the frequency
d = np.diff(r)
start = np.nonzero(d >= 0)[0][0]
peak_m = np.max(r[start:])
peak_i = np.argmax(r[start:])
peak_i = peak_i+start  # location of the 1st peak in the autocorrelation

plt.plot(peak_i, 0, 'rs', markersize=20)
plt.title ('Autocorrelation of spectrum')
plt.xlabel('lag, sample number')
plt.autoscale(enable=True, axis='x', tight=True)
plt.autoscale(enable=True, axis='y', tight=True)
plt.show

fsr = peak_i * np.mean(np.diff(lam1))
ng_av = np.mean(lam1)**2/(dL*1e-6)/fsr

print ('fsr =', fsr, 'ng_av = ', ng_av)

#==============================================================================
#  find starting point for curve fitting, using the ng value
#==============================================================================
lam0 = np.mean(lam1) * 1e6 # lambda0 is in microns.
n1=2.4
#modeNumber = n1_initial * dL / lambda0 - 0.5;
#n1 = (2*floor(modeNumber)+1)*lambda0/2/dL;
n2 = (n1-ng_av)/lam0
nx_init = np.array([n1, n2, 0])
alpha_init = 1e-3 # propagation loss [micron^-1]
x0 = np.hstack((nx_init, alpha_init, 0))

#==============================================================================
# Define the MZI transfer function
# - as a Taylor expansion around the central wavelength
# - Use units of [microns] - keeps the variables closer to 1.
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
print ('r2 = ', r2)

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
