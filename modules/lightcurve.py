"""
Module for analysing lightcurve data, specifically from TESS, although these functions are probably
    appropriate for lightcurves obtained by any telescope.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def plot_raw_flux(file):
    try:
        data = fits.open(file)
    except:
        raise Exception("failed to open fits file")
    header = data[1].header
    data = data[1].data

    data = data[~np.isnan(data['PDCSAP_FLUX'])]
    plt.scatter(data['TIME'], data['PDCSAP_FLUX']/np.mean(data['PDCSAP_FLUX']), s=0.05)
    plt.xlabel('Time (days)')
    plt.ylabel('Flux (e-/s)')
    plt.show()
