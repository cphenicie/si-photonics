# Python 2.7 script
# by Lukas Chrostowski in Matlab
# by Huixu (Chandler) Deng in Python, 2017

from __future__ import print_function # make 'print' compatible in Python 2X and 3X
from scipy.io import loadmat        # used to load MATLAB data files
import matplotlib.pyplot as plt     # used to plot data
import os                           # used to process filenames
import urllib                       # used to download files from Dropbox.com
import numpy as np

## Enter the order of the polynomial for the curve fitting:
#POLY_ORDER = 3

# Download the file from Dropbox.  Dropbox requires that you have a ?dl=1 at the end of the file
# Store the file in the local directory
# Enter the Dropbox URL here.  Make sure it has a =1 at the end:
#==============================================================================
# #  Loopback structure:
#==============================================================================
url_loopback = 'https://www.dropbox.com/s/w915qfix9kwlwv7/ZiheGao_MZI1_272_Scan1.mat?dl=1'
FileName_loopback = os.path.split(os.path.splitext(url_loopback)[0]+'.mat')[1]
print (FileName_loopback)

if not(os.path.exists(FileName_loopback)):
    print('Downloading from Dropbox')
    if sys.version_info[0] < 3:
        import urllib                       # used to download files from Dropbox.com
        urllib.urlretrieve (url, FileName_loopback) # used in Python 2.x
    else:
        from urllib.request import urlretrieve  # used to download files from Dropbox.com
        urlretrieve (url, FileName_loopback) # used in Python 3.x

PORT=1 # Which Fibre array port is the output connected to?
# Load the file, from the MATLAB format.
matData_loopback = loadmat(FileName_loopback, squeeze_me=True, struct_as_record=False)
# Read the experimental data from the MATLAB file
lam_loopback = matData_loopback['scanResults'][PORT-1].Data[:,0]/1e9
amplitude_loopback = matData_loopback['scanResults'][PORT-1].Data[:,1]

#==============================================================================
# #  MZI:
#==============================================================================
url_mzi = 'https://www.dropbox.com/s/1rvjfef4jqybc12/ZiheGao_MZI2_271_Scan1.mat?dl=1'
FileName_mzi = os.path.split(os.path.splitext(url_mzi)[0]+'.mat')[1]

if not(os.path.exists(FileName_mzi)):
    print('Downloading from Dropbox')
    if sys.version_info[0] < 3:
        import urllib                       # used to download files from Dropbox.com
        urllib.urlretrieve (url, FileName_mzi) # used in Python 2.x
    else:
        from urllib.request import urlretrieve  # used to download files from Dropbox.com
        urlretrieve (url, FileName_mzi) # used in Python 3.x

PORT=1 # Which Fibre array port is the output connected to?
# Load the file, from the MATLAB format.
matData_mzi = loadmat(FileName_mzi, squeeze_me=True, struct_as_record=False)
# Read the experimental data from the MATLAB file
lam_mzi = matData_mzi['scanResults'][PORT-1].Data[:,0]/1e9
amplitude_mzi = matData_mzi['scanResults'][PORT-1].Data[:,1]

#==============================================================================
# # Calibrate the MZI data using the loopback structure
# # Plot
#==============================================================================
FONTSIZE=20

plt.figure()
plt.plot (lam_loopback*1e6, amplitude_loopback)
plt.title ('Calibration loopback',fontsize = FONTSIZE)
plt.xlabel ('Wavelength [$\mu$m]',fontsize = FONTSIZE)
plt.ylabel ('Insertion Loss [dB]',fontsize = FONTSIZE)

# Fit the data with a polynomial
p_loopback=np.polyfit((lam_loopback-np.mean(lam_loopback))*1e6, amplitude_loopback, 5)
amplitude_loopback_poly=np.polyval(p_loopback,(lam_loopback-np.mean(lam_loopback))*1e6)
plt.plot(lam_loopback*1e6, amplitude_loopback_poly, 'r')

# find wavelength range with usable data, in the loopback
loopback_IL = np.max(amplitude_loopback)
new_lam_loopback_i=np.nonzero(amplitude_loopback > loopback_IL-10)
lam_loopback=lam_loopback[new_lam_loopback_i]
lam_loopback_min = np.min(lam_loopback)
lam_loopback_max = np.max(lam_loopback)
amplitude_loopback=amplitude_loopback[new_lam_loopback_i]

# refit the loopback
LOOPBACK = np.polyfit((lam_loopback-np.mean(lam_loopback))*1e6, amplitude_loopback, 4)
amplitude_loopback_poly = np.polyval(LOOPBACK,(lam_loopback-np.mean(lam_loopback))*1e6)

plt.plot(lam_loopback*1e6, amplitude_loopback_poly,'r-', linewidth = 5)
plt.autoscale(enable=True, axis='x', tight=True)
plt.autoscale(enable=True, axis='y', tight=True)
plt.show()

#==============================================================================
# MZI data:
#==============================================================================
plt.figure()
plt.plot (lam_mzi*1e6, amplitude_mzi)
plt.title ('MZI (raw data)',fontsize = FONTSIZE)
plt.xlabel ('Wavelength [$\mu$m]',fontsize = FONTSIZE)
plt.ylabel ('Insertion Loss [dB]',fontsize = FONTSIZE)
plt.show()

#==============================================================================
# MZI data - calibrated
#==============================================================================
# data only within the bandwidth of interest.
lam_mzi_p=np.arange(lam_loopback_min, lam_loopback_max, np.min(np.diff(lam_mzi)))
amplitude_mzi_p=np.interp(lam_mzi_p, lam_mzi, amplitude_mzi)
amplitude_mzi_p[np.nonzero(amplitude_mzi_p==-np.inf)]=-50
                
# calibrate data
amplitude_mzi_cal=amplitude_mzi_p-np.polyval(LOOPBACK,(lam_mzi_p-np.mean(lam_mzi_p))*1e6);
plt.figure()
plt.plot(lam_mzi_p*1e6, amplitude_mzi_cal);
plt.title ('MZI (calibrated with loopback',fontsize = FONTSIZE)
plt.xlabel ('Wavelength [$\mu$m]',fontsize = FONTSIZE)
plt.ylabel ('Insertion Loss [dB]',fontsize = FONTSIZE)
plt.show()
