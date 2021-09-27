"""
Plot raw values of flux vs time from a directory of fits files containing
    lightcurve data for a given system.
Useful when scraping a large number of lightcurves from the archive.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

ls = os.listdir()
files = [('.fits' in file) & ('llc' in file) for file in ls]
#files = ['.fits' in file for file in ls]
files = np.asarray(ls)[files]

for file in files:
    try:
        data = fits.open(file)
    except:
        continue
    header = data[1].header
    data = data[1].data

    data = data[~np.isnan(data['PDCSAP_FLUX'])]

    plt.plot(data['TIME'], data['PDCSAP_FLUX']/np.mean(data['PDCSAP_FLUX']))#, s=0.01)
    plt.xlabel('Time (days)')
    plt.ylabel('Flux (e-/s)')
plt.show()
