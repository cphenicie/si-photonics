# Python 2.7 script
# by Lukas Chrostowski in Matlab
# by Huixu (Chandler) Deng in Python, 2017

from __future__ import print_function # make 'print' compatible in Python 2X and 3X
from scipy.io import loadmat        # used to load MATLAB data files
import matplotlib.pyplot as plt     # used to plot data
import os                           # used to process filenames
import sys

# Download the file from Dropbox.  Dropbox requires that you have a ?dl=1 at the end of the file
# Store the file in the local directory
# Enter the Dropbox URL here.  Make sure it has a =1 at the end:
url = 'https://www.dropbox.com/s/1rvjfef4jqybc12/ZiheGao_MZI2_272_Scan1.mat?dl=1'
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
wavelength = matData['scanResults'][PORT-1].Data[:,0]/1e9
power = matData['scanResults'][PORT-1].Data[:,1]

# Plot the raw data:
plt.figure()
plt.plot(wavelength*1e6, power)
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('Transmission [dB]')
plt.autoscale(enable=True, axis='x', tight=True)
plt.autoscale(enable=True, axis='y', tight=True)
plt.title('Experimental data (raw)')