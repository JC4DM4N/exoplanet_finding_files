"""
Generate output files from fits files
Output file has columns:
    time, noramlised flux, normalised error
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

ls = os.listdir()
files = [('.fits' in file) & ('llc' in file) for file in ls]
#files = ['.fits' in file for file in ls]
files = np.asarray(ls)[files]

ndays = []
ifile = []

for i, file in enumerate(files):
    try:
        data = fits.open(file)
    except:
        continue
    header = data[1].header
    data = data[1].data

    data = data[~np.isnan(data['PDCSAP_FLUX'])]
    time = data['TIME']
    # convert to Julian days
    time += 2454833
    flux = data['PDCSAP_FLUX']
    fluxerr = data['PDCSAP_FLUX_ERR']

    fluxnorm = flux/np.mean(flux)
    fluxerrnorm = fluxerr/flux*fluxnorm

    outfile = open('output%i.dat' %i, 'w')
    for j in range(len(time)):
        if j==0:
            t0 = time[j]
        outfile.write('%s %s %s 0 0 -1 \n' %(str(time[j]), str(fluxnorm[j]), str(fluxerrnorm[j])))
    ifile.append(i)
    ndays.append(time[j] - t0)
    outfile.close()

print('Number of days of data in each file: ')
print(ifile)
print(ndays)
