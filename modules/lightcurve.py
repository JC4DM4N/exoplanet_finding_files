"""
Module for analysing lightcurve data, specifically from TESS, although these functions are probably
    appropriate for lightcurves obtained by any telescope.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def read_flux_vs_time(file,flux_label='SAP_FLUX'):
    """
    Read fits file and return time & flux values
    Inputs:
        flux_label : may either be SAP_FLUX or PDCSAP_FLUX, depending whether GP removed
            or GP included data is wanted.
    """
    try:
        data = fits.open(file)
    except:
        raise Exception("failed to open fits file")
    header = data[1].header
    data = data[1].data
    time, flux = remove_nans(data['TIME'], data[flux_label])
    return time, flux

def remove_nans(time, flux):
    """
    remove data points where flux is NaN.
    """
    mask = ~np.isnan(flux)
    return time[mask], flux[mask]

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
    time_BJD = time + offset
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

def bin_data(time,flux,bin_size=30):
    """
    Function to bin data every X minutes
    Inputs:
        bin_size : duration of bins in minutes
    """
    # convert bin size to days
    bin_size = bin_size/60./24.
    time_bins = np.arange(np.min(time), np.max(time), bin_size)
    ibins = np.digitize(time, time_bins)
    binned_fluxes = np.asarray([np.mean(flux[ibins==i]) for i in np.arange(ibins.max())])
    return time_bins, binned_fluxes

def plot_raw_flux(file):
    """
    Open fits file and plot PDCSAP_FLUX vs TIME.
    Inputs:
        file : path to fits file.
    """
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
