"""
Read the RV file and plot values of RV vs time with associated RV errors.
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
sys.path.insert(0,os.path.join(os.path.dirname(os.path.realpath(__file__)),'../../modules'))
import radialvelocity as rv

parser = argparse.ArgumentParser()
parser.add_argument('-file')
args = parser.parse_args()
file = args.file

# load the data file with the RVs and their errors and plot the RVs vs time
rv.plot_raw_rvs(file, save=True, show=True)
plt.show()
