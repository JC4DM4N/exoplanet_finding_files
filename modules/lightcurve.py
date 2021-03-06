"""
Module for analysing lightcurve data, specifically from TESS, although these functions are probably
    appropriate for lightcurves obtained by any telescope.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
from astropy.io import fits

def read_flux_vs_time(file,flux_label='PDCSAP_FLUX'):
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
    time, flux, err = remove_nans(data['TIME'], data[flux_label], data[flux_label+"_ERR"])
    return time, flux, err

def remove_nans(time, flux, flux_err):
    """
    remove data points where flux is NaN.
    """
    mask = ~np.isnan(flux)
    return time[mask], flux[mask], flux_err[mask]

def phasefold_data(times,period=6.51595066707795):
    """
    Fold the time values on the planet transit period.
    """
    return np.mod(times,period)

def remove_trasit_points(time,flux,err,period=6.51595066707795,duration=2.53343092705777,offset=2457000,Tc=2458876.00,plot=True):
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
    return time[~mask], flux[~mask], err[~mask]

def bin_data(time,flux,err,bin_size=30):
    """
    Function to bin data every X minutes
    Inputs:
        bin_size : duration of bins in minutes
    """
    # convert bin size to days
    bin_size = bin_size/60./24.
    time_bins = np.arange(np.min(time), np.max(time), bin_size)
    ibins = np.digitize(time, time_bins)
    binned_fluxes = []
    binned_flux_errs = []
    for bin in np.unique(ibins):
        mask = ibins==bin
        binned_flux = np.mean(flux[mask])
        binned_fluxes.append(binned_flux)
        binned_flux_err = binned_flux*np.sqrt(np.sum(np.square(err[mask]/flux[mask])))
        binned_flux_err = binned_flux_err/np.sum(mask)
        binned_flux_errs.append(binned_flux_err)
    time_bins = time_bins[np.unique(ibins)-1]
    binned_fluxes = np.asarray(binned_fluxes)
    binned_flux_errs = np.asarray(binned_flux_errs)

    # remove nans, which correpsond to empty  bins
    return remove_nans(time_bins, binned_fluxes, binned_flux_errs)

def convert_to_BJD(time,offset=2457000):
    return time+offset

def subtract_mean_flux(flux):
    """
    Return mean subtracted fluxes from raw fluxes.
    """
    flux -= np.median(flux)
    return flux

def plot_raw_flux(file,flux_label='PDCSAP_FLUX',save=False):
    """
    Open fits file and plot PDCSAP_FLUX vs TIME.
    Inputs:
        file : path to fits file.
    """
    # read data
    t, f, e = read_flux_vs_time(file,flux_label)
    plt.figure(figsize=(10,5))
    plt.scatter(t, f/np.mean(f), s=1, c='black')
    plt.xlabel('Time (BJD - 2457000 days)',fontsize=15)
    plt.ylabel('Relative Flux (e-/s)',fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.4f}'))
    plt.xlim([min(t),max(t)])
    if save:
        plt.tight_layout()
        plt.savefig('TESS_LC_flux_vs_time.pdf')
    plt.show()

def plot_binned_flux(file,flux_label='PDCSAP_FLUX',save=False):
    """
    Open fits file and plot PDCSAP_FLUX vs TIME.
    Inputs:
        file : path to fits file.
    """
    # read data
    t, f, e = read_flux_vs_time(file,flux_label)
    # bin data
    tbin, fbin, ebin = bin_data(t,f,e)
    plt.figure(figsize=(10,5))
    plt.scatter(tbin, fbin/np.mean(fbin), s=5, c='black')
    plt.xlabel('Time (BJD - 2457000 days)',fontsize=15)
    plt.ylabel('Relative Flux (e-/s)',fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.4f}'))
    plt.xlim([min(tbin),max(tbin)])
    if save:
        plt.tight_layout()
        plt.savefig('TESS_LC_flux_vs_time_%s.pdf' %flux_label)
    plt.show()

def plot_binned_phase_folded_flux(file,period=6.51595066707795,Tc=2458876.025023,duration=2.53343092705777+0.6371884,save=False):
    """
    Open fits file and plot PDCSAP_FLUX vs TIME.
    Inputs:
        file : path to fits file.
    """
    # read data
    t, f, e = read_flux_vs_time(file,'PDCSAP_FLUX')
    # phasefold data
    times_pf = phasefold_data(t,period)
    # bin data
    tbin, fbin, ebin = bin_data(times_pf,f,e)

    # now zoom in on region of the transit centre
    # get relative position of tranist centre, and plot a few hours either side of this
    RTc = np.mod(Tc-2457000,period)
    transit = (tbin > RTc - duration/2./24.) & (tbin < RTc + duration/2./24)
    xmin, xmax = (RTc - period/10./2., RTc + period/10./2.) # plot 1/10 of a phase
    # convert times from days to phase
    tbin = tbin/period
    times_pf = times_pf/period
    xmin = xmin/period
    xmax = xmax/period

    plt.figure(figsize=(10,5))
    plt.errorbar(times_pf,f/np.mean(f),yerr=e/np.mean(f),fmt='.',c='grey',zorder=0)
    plt.errorbar(tbin[~transit],fbin[~transit]/np.mean(fbin),yerr=ebin[~transit]/np.mean(fbin),fmt='.',c='black',zorder=1)
    plt.errorbar(tbin[transit],fbin[transit]/np.mean(fbin),yerr=ebin[transit]/np.mean(fbin),fmt='.',c='red',zorder=2)
    plt.xlabel('Phase',fontsize=15)
    plt.ylabel('Relative Flux (e-/s)',fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlim([0,1])
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.4f}'))
    if save:
        plt.tight_layout()
        plt.savefig('TESS_LC_pf_flux_vs_time.pdf')

    plt.figure(figsize=(10,5))
    plt.errorbar(times_pf,f/np.mean(f),yerr=e/np.mean(f),fmt='.',c='grey',zorder=0)
    plt.errorbar(tbin[~transit],fbin[~transit]/np.mean(fbin),yerr=ebin[~transit]/np.mean(fbin),fmt='.',c='black',zorder=1)
    plt.errorbar(tbin[transit],fbin[transit]/np.mean(fbin),yerr=ebin[transit]/np.mean(fbin),fmt='.',c='red',zorder=1)
    plt.xlabel('Phase',fontsize=15)
    plt.ylabel('Relative Flux (e-/s)',fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlim([xmin,xmax])
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.4f}'))
    if save:
        plt.tight_layout()
        plt.savefig('TESS_LC_pf_zoomed_flux_vs_time.pdf')
    plt.show()
