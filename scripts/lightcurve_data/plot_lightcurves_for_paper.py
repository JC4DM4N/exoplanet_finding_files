"""
Generate plots of flux vs time for both PDCSAP_FLUX and SAP_FLUX data.
Fluxes are binned every 30 mins.
"""

import sys
import os
import argparse
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0,os.path.join(os.path.dirname(os.path.realpath(__file__)),'../../modules'))
import lightcurve as lc

parser = argparse.ArgumentParser()
parser.add_argument('-file') # path to raw fits file
args = parser.parse_args()
file = args.file

for label in ['PDCSAP_FLUX','SAP_FLUX']:
    lc.plot_binned_flux(file,flux_label=label,save=True)
lc.plot_phase_folded_flux(file,save=True)
