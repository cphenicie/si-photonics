# Python 3 script
#
# - Download a MATLAB measurements data file from the web (Dropbox)
# - Plot the data, save to a PDF file
#
# by Lukas Chrostowski, 2015
#

from scipy.io import loadmat      # used to load MATLAB data files
import matplotlib.pyplot as plt   # used to plot data import os # used to process filenames
import urllib                     # used to download files from Dropbox.com
import os

# Download the file from Dropbox. Dropbox requires that you have a ?dl=1 at the end of the file
# Store the file in the local directory
url = "https://www.dropbox.com/s/1rvjfef4jqybc12/ZiheGao_MZI2_271_Scan1.mat?dl=1"
FileName = os.path.split(os.path.splitext(url)[0]+'.mat')[1]
print (FileName)
# urllib.urlretrieve (url, FileName)
urllib.request.urlretrieve(url, FileName)

# Load the file, from the MATLAB format.
matData = loadmat(FileName, squeeze_me=True, struct_as_record=False)

# Read the experimental data from the MATLAB file
wavelength = matData['scanResults'][0].Data[:,0]/1e9
power = matData['scanResults'][0].Data[:,1]

# Create a figure; plot the data; save the figure to a PDF file
plt.figure()
plt.plot(wavelength*1e9, power)
plt.xlim((1530,1560))
plt.xlabel('Wavelength (nm)')
plt.ylabel('Transmission (dB)')
saveFigFileName = os.path.splitext(FileName)[0]+'.pdf'
plt.title('Raw spectrum data: %s'%(saveFigFileName))
plt.savefig(saveFigFileName)
saveFigFileName = os.path.splitext(FileName)[0]+'.png'
plt.title('Raw spectrum data: %s'%(saveFigFileName))
plt.savefig(saveFigFileName)
plt.show()