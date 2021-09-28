"""
Module for analysing lightcurve data, specifically from TESS, although these functions are probably
    appropriate for lightcurves obtained by any telescope.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def read_flux_vs_time(file):
    """
    Read fits file and return time & flux values
    """
    try:
        data = fits.open(file)
    except:
        raise Exception("failed to open fits file")
    header = data[1].header
    data = data[1].data
    return data['TIME'], data['PDCSAP_FLUX']

def phasefold_data(times,period):
    """
    Fold the time values on the planet transit period.
    """
    return np.mod(times,period)

def remove_trasit_points(time,flux,period=6.51595066707795,duration=2.53343092705777,offset=-2458870.4338364583):
    """
    Inputs:
        period : period of the transiting planet (days).
        duration : duration of the transit (hours).
        offset : represents the time since time[0] that the first transit begins, such that transits
            occur at time[0] + offset + N*period (days).
    Outputs:
        time, flux : with points removed which are during a transit.
    """
    # first phase fold the data
    time = phasefold_data(time,period)

    #mask = (time >= time[0] + offset + period) & (time <= time[0] + offset + period + duration)
    import pdb; pdb.set_trace()
    return


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
