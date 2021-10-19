import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-file')
args = parser.parse_args()
file = args.file

data = np.genfromtxt(file, delimiter=' ')
#want roughly 3 bins per day
nbins = int((data[-1,0] - data[0,0])*3.)
timebins = np.histogram(data[:,0],nbins)[1]
ibins = np.digitize(data[:,0],bins=timebins)
ibins -= 1 # change so index starts at 0

outfile = open('%s_binned.dat' %(file.split('.dat')[0]), 'w')

for i, bin in enumerate(timebins):
    #want the mean flux and the error on the mean for each time bin
    inbin = ibins == i

    if np.sum(inbin) == 0:
        continue

    #inbin = (data[:,0] >= timebins[i]) & (data[:,0] < timebins[i])
    fluxes = data[inbin,1]
    err = data[inbin,2]
    errsum = 0
    #sum the errors in quadrature
    for j, flux in enumerate(fluxes):
        errsum += (err[j])**2

    timemean = np.mean(data[inbin,0])
    fluxmean = np.mean(fluxes)
    #error on the mean
    errmean = np.sqrt(errsum)/len(err)

    outfile.write('%s %s %s 0 0 -1 \n' %(timemean, fluxmean, errmean))

outfile.close()
