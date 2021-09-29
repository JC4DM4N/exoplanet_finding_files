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

def remove_trasit_points(time,flux,period=6.51595066707795,duration=2.53343092705777,offset=2457000,Tc=2458876.00,plot=True):
    """
    Inputs:
        period : period of the transiting planet (days).
        duration : duration of the transit (hours).
        offset : time array offset from BJD=0.
        Tc : represents the time of transit centre, required to determine where the transit
            is in the lightcurve.
    Outputs:
        time, flux : with points removed which are during a transit.
    """
    # convert times to BJD
    time_BJD += offset
    # first phase fold the data
    time_pf = phasefold_data(time_BJD,period)
    Tc_pf = phasefold_data(Tc,period)
    # remove points during the transit
    duration = duration/24.
    mask = (time_pf >= Tc_pf - duration/2.) & (time_pf <= Tc_pf + duration/2.)

    if plot:
        plt.figure('original data')
        plt.scatter(time_pf,flux,s=0.1)
        plt.figure('transit removed')
        plt.scatter(time_pf[~mask],flux[~mask],s=0.1)
        plt.show()
    return time[~mask], flux[~mask]


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
