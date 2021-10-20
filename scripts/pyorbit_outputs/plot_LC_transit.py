"""
Plot best fit transit along with raw flux data and binned flux data for tranist fit
    output file from PyORBIT.
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

#LC_data_full = np.genfromtxt('LCdata_full.dat',skip_header=1)
data = np.genfromtxt('LCdata_lc_model_b.dat',skip_header=1)

# lets bin the LC data every 30 mins
# assume phase is 6.52 days, hence 30 mins is phase = time/6.52days
time_bins = np.arange(np.min(data[:,2]),np.max(data[:,0]),30./60./24./6.52)
ibins = np.digitize(data[:,2], time_bins)
phase_binned = []
flux_binned = []
flux_err_binned = []
residuals_binned = []
residuals_err_binned = []
for bin in np.unique(ibins):
    mask = ibins==bin
    phase_binned.append(np.mean(data[mask,2]))
    flux = np.mean(data[mask,3])
    # sum all the relative errors in quadrature, then square root and multiply by the mean flux
    err = flux*np.sqrt(np.sum(np.square(data[mask,4]/data[mask,3])))
    # error on the mean - take a factor 1/N
    err = err/(np.sum(mask))
    flux_binned.append(flux)
    flux_err_binned.append(err)

    # also bin the residuals and their errors
    residuals = np.mean(data[mask,10])
    residuals_err = residuals*np.sqrt(np.sum(np.square(data[mask,11]/data[mask,10])))
    residuals_err = residuals_err/(np.sum(mask))
    residuals_binned.append(residuals)
    residuals_err_binned.append(residuals_err)

plt.figure(figsize=(6,5))

gs = gridspec.GridSpec(2,1,height_ratios=[3,1])
plt.subplot(gs[0])
# scatter plot the raw flux values
plt.errorbar(data[:,2],data[:,3],yerr=data[:,4],fmt='.',color='gray',zorder=0,elinewidth=0.5)
# scatter plot the best fit transit fluxes
plt.scatter(data[:,2],data[:,7],color='red',s=0.5,zorder=1)
# also scatter plot the binned flux data
plt.errorbar(phase_binned,flux_binned,yerr=flux_err_binned,c='black',fmt='.',zorder=2,elinewidth=0.5)
plt.ylabel('Relative flux (e-/s)')
plt.xlim(0.164,0.27)

plt.subplot(gs[1])
# now plot the residuals
plt.errorbar(data[:,2],data[:,10],yerr=data[:,11],fmt='.',color='gray',elinewidth=0.5,zorder=0)
# and overplot the binned residuals
plt.errorbar(phase_binned,residuals_binned,yerr=residuals_err_binned,fmt='.',color='black',elinewidth=0.5,zorder=1)
plt.xlabel('Phase')
plt.ylabel('Residuals')
plt.xlim(0.164,0.27)

plt.tight_layout()

plt.savefig('LC_transit_TOI1778.png')
plt.show()
