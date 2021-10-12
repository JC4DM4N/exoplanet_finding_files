"""
Read the RV file and plot values of RV vs time with associated RV errors.
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-file')
args = parser.parse_args()
file = args.file

# load the data file with the RVs and their errors
data = np.genfromtxt(file)
time = data[:,0]
rv = data[:,1]
rv_error = data[:,2]

#plot the RVs vs time
plt.errorbar(time,rv,yerr=rv_error,fmt='.',c='black',linewidth=1)
plt.ylabel(r'RV (ms$^{-1}$)',fontsize=15)
plt.xlabel('BJD - 2450000 (Days)',fontsize=15)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.tight_layout()
plt.savefig('raw_RV_plot.png')
plt.show()
