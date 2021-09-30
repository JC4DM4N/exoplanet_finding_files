"""
Prepare raw lightcurve data for input to PyORBIT GP analysis.
Take TESS lightcurve which contains a transit, remove all data points which
    represent a transit, and bin the remaining data every 30 mins.
"""

import sys
import os
import argparse
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0,os.path.join(os.path.dirname(os.path.realpath(__file__)),'../modules'))
import lightcurve as lc

parser = argparse.ArgumentParser()
parser.add_argument('-file') # path to raw fits file
parser.add_argument('-plot',action='store_true') # boolean, where to plot results
args = parser.parse_args()
file = args.file
plot = args.plot

# load data
t,f = lc.read_flux_vs_time(file)
# remove data at transit
t,f = lc.remove_trasit_points(t,f,plot=False)
# bin the remaining data every 30 mins
tbin,fbin = lc.bin_data(t,f)

if plot:
    plt.scatter(tbin,fbin,s=0.1)
    plt.xlabel('Time')
    plt.ylabel('Flux')
    plt.show()

with open('TESS_LC.dat','w') as file:
    for t,f in zip(tbin,fbin):
        file.write('%.8f %.3f \n' %(t,f))
