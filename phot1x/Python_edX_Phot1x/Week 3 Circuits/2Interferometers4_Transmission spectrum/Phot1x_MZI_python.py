# Python 2.7 script
# by Lukas Chrostowski in Matlab
# by Huixu (Chandler) Deng in Python, 2017

from __future__ import print_function # make 'print' compatible in Python 2X and 3X
import matplotlib.pyplot as plt
import numpy as np

# the wavelength range of interest.
lam_min = 1.5  # Units [µm, microns]
lam_max = 1.6
lam_step = 0.01e-3 # wavelength step [microns], Typical minimum step for a tunable laser is 1-10 pm.
lam=np.arange(lam_min, lam_max, lam_step)

# Effective index:
# - as a Taylor expansion around the central wavelength, lam0
lam0 = 1.55
n1=2.4
n2=-1.0
n3=0.0  # these are constants from the waveguide model.

def neff(lam):
    return n1 + n2*(lam-lam0) + n3*(lam-lam0)**2

# plot, and check if this is as expected:
plt.figure()
plt.plot(lam, neff(lam),linewidth=3.0)

# Complex propagation constant
alpha = 1e-3   # propagation loss [micron^-1]; constant

def beta(lam):
    return	2*np.pi*neff(lam)/lam - 1j*alpha/2*np.ones(np.size(lam))
    
# MZI transfer function
def T_MZI(L1, L2, lam):
    return 0.25*np.abs(np.exp(-1j*beta(lam)*L1) + np.exp(-1j*beta(lam)*L2))**2

# plot, and check if this is as expected:
L1=100
L2=200  # Units [µm, microns], variable

plt.figure()
plt.plot(lam, T_MZI(L1, L2, lam), linewidth=3)
plt.xlabel ('Wavelength [$\mu$m]')
plt.ylabel ('Transmission')
plt.autoscale(enable=True, axis='x', tight=True)
plt.autoscale(enable=True, axis='y', tight=True)
plt.title('MZI transfer function')

plt.figure()
T_MZI_dB = 10*np.log10(T_MZI(L1, L2, lam))
plt.plot(lam, T_MZI_dB, linewidth=3)
plt.xlabel ('Wavelength [$\mu$m]')
plt.ylabel ('Transmission [dB]')
plt.autoscale(enable=True, axis='x', tight=True)
plt.autoscale(enable=True, axis='y', tight=True)
plt.title ('MZI transfer function')
