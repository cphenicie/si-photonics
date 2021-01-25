from siepicfab import all as pdk
from ipkiss3 import all as i3
import numpy as np
import pylab as plt

splitter = pdk.EbeamY1550()

# 1. Layout
splitter_layout = splitter.Layout()
splitter_layout.visualize(annotate=True)
splitter_layout.visualize_2d()
splitter_layout.cross_section(cross_section_path=i3.Shape([(-0.5, -1.5), (-0.5, 1.5)])).visualize()

# 2. Circuit
splitter_circuit = splitter.CircuitModel()

# 3. Plotting
wavelengths = np.linspace(1.5, 1.58, 51)
S = splitter_circuit.get_smatrix(wavelengths=wavelengths)
plt.plot(wavelengths, 10 * np.log10(np.abs(S["opt2", "opt1"]) ** 2), '-', linewidth=2.2, label="T(opt2)")
plt.plot(wavelengths, 10 * np.log10(np.abs(S["opt3", "opt1"]) ** 2), '-', linewidth=2.2, label="T(opt1)")
plt.ylim(-5., 0.)
plt.xlabel('Wavelength [um]', fontsize=16)
plt.ylabel('Transmission [dB]', fontsize=16)
plt.legend(fontsize=14, loc=1)
plt.show()

