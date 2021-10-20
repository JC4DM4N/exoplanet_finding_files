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
    omegas = 2. * pi * f
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
