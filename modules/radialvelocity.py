"""
Module for analysing radial velocity data, specifically from HARPS-N, although these functions are probably
    appropriate for RVs obtained by any spectrograph.
"""

import numpy as np
import matplotlib.pyplot as plt

def bgls(t, y, err, plow=1.0, phigh=100.0, ofac=10, jit=0.0, dt = None):
    """
    BGLS: Calculates Bayesian General Lomb-Scargle periodogram, normalised with minimum.

    t: times
    y: data
    err: error bars
    plow: lowest period to sample
    phigh: highest period to sample
    ofac: oversampling factor
    jit: white noise to be added to the error bars
    dt: time span to be considered
    """
    # Define time span
    if dt == None:
        dt = np.max(t)-np.min(t)

    # Amount of frequencies to sample
    amount = (phigh-plow)*dt*ofac/plow/phigh + 1
    # Define frequencies
    f = np.linspace(1./phigh, 1./plow, int(amount))
    omegas = 2.*np.pi*f
    # Define weights - optional white noise can be added
    err2 = err * err + jit * jit
    w = 1./err2
    W = sum(w)
    # Follow algorithm in Mortier et al. 2015
    bigY = sum(w*y)

    p = []
    constants = []
    exponents = []
    for i, omega in enumerate(omegas):
        theta = 0.5 * np.arctan2(sum(w*np.sin(2.*omega*t)), sum(w*np.cos(2.*omega*t)))
        x = omega*t - theta
        cosx = np.cos(x)
        sinx = np.sin(x)
        wcosx = w*cosx
        wsinx = w*sinx
        C = sum(wcosx)
        S = sum(wsinx)
        YCh = sum(y*wcosx)
        YSh = sum(y*wsinx)
        CCh = sum(wcosx*cosx)
        SSh = sum(wsinx*sinx)
        if (CCh != 0 and SSh != 0):
            K = ((C*C*SSh + S*S*CCh - W*CCh*SSh)/(2.*CCh*SSh))
            L = ((bigY*CCh*SSh - C*YCh*SSh - S*YSh*CCh)/(CCh*SSh))
            M = ((YCh*YCh*SSh + YSh*YSh*CCh)/(2.*CCh*SSh))
            constants.append(1./np.sqrt(CCh*SSh*abs(K)))
        elif (CCh == 0):
            K = (S*S - W*SSh)/(2.*SSh)
            L = (bigY*SSh - S*YSh)/(SSh)
            M = (YSh*YSh)/(2.*SSh)
            constants.append(1./np.sqrt(SSh*abs(K)))
        elif (SSh == 0):
            K = (C*C - W*CCh)/(2.*CCh)
            L = (bigY*CCh - C*YCh)/(CCh)
            M = (YCh*YCh)/(2.*CCh)
            constants.append(1./np.sqrt(CCh*abs(K)))
        if K > 0:
            raise RuntimeError('K is positive. This should not happen.')
        exponents.append(M - L*L/(4.*K))

    constants = np.array(constants)
    exponents = np.array(exponents)
    logp = (np.log10(constants) + (exponents * np.log10(np.exp(1.))))
    # Normalise
    logp = logp - min(logp)
    # Return array of frequencies and log of probability
    return f, logp

def load_rvs(file):
    """
    Load RVs from a datafile with columns [TIME, RV, RV_ERR]
    This is the same file format as is used as inputs for the PyORBIT analysis.
    """
    # load the data file with the RVs and their errors
    data = np.genfromtxt(file)
    time = data[:,0]
    rv = data[:,1]
    rv_error = data[:,2]
    return time, rv, rv_error

def plot_raw_rvs(file, show=False, save=False):
    """
    Load RVs from a data file, and plot time vs. RV with associated error bars.
    """
    time, rv, rv_err = load_rvs(file)

    #plot the RVs vs time
    plt.errorbar(time,rv,yerr=rv_error,fmt='.',c='black',linewidth=1)
    plt.ylabel(r'RV (ms$^{-1}$)',fontsize=15)
    plt.xlabel('BJD - 2450000 (Days)',fontsize=15)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    plt.tight_layout()
    if save:
        plt.savefig('raw_RV_plot.png')
    if show:
        plt.show()
