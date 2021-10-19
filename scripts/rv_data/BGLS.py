# Copyright (c) 2016-2017 A. Mortier
# Distributed under the MIT License
###################################
# Permission is hereby granted, free of charge,
# to any person obtaining a copy of this software
# and associated documentation files (the "Software"),
# to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify,
# merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to
# whom the Software is furnished to do so, subject
# to the following conditions:

# The above copyright notice and this permission
# notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT
# WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
# ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
##################################

# This code was the basis for the following works
# Mortier et al. 2015
# http://adsabs.harvard.edu/abs/2015A%26A...573A.101M
# Mortier et al. 2017
# http://adsabs.harvard.edu/abs/2017A%26A...601A.110M
# Please cite these works if you use this.

# Import necessary packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--file','-file',help='path to rdb file')
parser.add_argument('--ydata','-ydata',help='column header for the y data (vrad/fwhm)')
args = parser.parse_args()
fn = args.file
ydata = args.ydata

# want to ensure this is the rdb file, which contains fwhm columns
assert fn.endswith('.rdb')
# rdb file has rows [header, placeholder, data], so we want to drop the placeholder row
data = pd.read_csv(fn,delimiter='\t')
data = data.drop(index=0)

# Predefine pi
pi = np.pi

t = np.asarray(data['rjd']).astype(float)
y = np.asarray(data[ydata]).astype(float)
if ydata=='vrad':
    err = np.asarray(data['svrad']).astype(float)
elif ydata=='fwhm':
    err = np.asarray(data['sig_fwhm']).astype(float)
else:
    raise Exception('Unfamiliar ydata header used. Bailing out.')

print(t, y, err)

def bgls(t, y, err, plow=1.0, phigh=100.0, ofac=10, jit=0.0, dt = None):
    """
    BGLS: Calculates Bayesian General Lomb-Scargle
                      periodogram, normalised with minimum.

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

f, logp = bgls(t,y,err)
print(f, logp)

plt.figure(figsize=(15,5))
plt.semilogx(1.0/f,logp,label='HARPSN : TOI-1778', linewidth=1.0)
plt.xlim(1.0,100.0)
plt.plot([6.526,6.526],[0.0,100.0],color='black',ls='dashed',lw=0.75)
plt.xlabel('Period (days)')
plt.ylabel('log$_{10}$ (Power)')
plt.legend(loc='upper center')
plt.ylim(0.0,100.0)
plt.savefig('TOI-1778_BGLS.png',dpi=1000)
plt.show()

def sbgls(t, y, err, obsstart = 5, plow=0.5, phigh=100,
          ofac=1, jit=0.0, fig='yes'):

    '''
    SBGLS: Calculates Stacked BGLS periodogram
           Can plot figure

    t: times
    y: data
    err: error bars
    obsstart: minimum number of observations to start
    plow: lowest period to sample
    phigh: highest period to sample
    ofac: oversampling factor
    jit: white noise to be added to the error bars
    fig: creating figure
    '''

    n = len(t)

    # Timespan
    dt = np.max(t)-np.min(t)

    # Empty lists to fill for sbgls
    freqs = []
    powers = []
    nrobs = []

    # Do BGLS for every set of observations
    # and save results
    for i in range(obsstart,n+1):

        freq, power = bgls(t[:i], y[:i], err[:i],
                           plow=plow, phigh=phigh,
                           ofac=ofac, jit=jit, dt = dt)

        freqs.extend(freq)
        powers.extend(power)
        nrobs.extend(np.zeros(len(freq))+i)

    freqs = np.array(freqs)
    powers = np.array(powers)
    nrobs = np.array(nrobs)

    # Make figure
    if fig == 'yes':

        plt.scatter(1./freqs, nrobs, c=powers,
                    cmap = plt.get_cmap('Reds'),
                    lw=0, marker = 's')

        plt.xlim(min(1./freqs),max(1./freqs))
        plt.ylim(min(nrobs),max(nrobs))

        plt.xlabel('Period (days)',fontsize=16)
        plt.ylabel('Nr of observations',fontsize=16)

        plt.gca().set_xscale('log')

        cbar = plt.colorbar()
        cbar.set_label(r'$\log P$',fontsize=14)

    return freqs, nrobs, powers

#freq, nrobs, powers = sbgls(t,y,err)
