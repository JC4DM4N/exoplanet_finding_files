"""
Plot the BGLS power spectrum for one or more variables.
"""

import sys
import os
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
sys.path.insert(0,os.path.join(os.path.dirname(os.path.realpath(__file__)),'../../modules'))
import radialvelocity as rv

parser = argparse.ArgumentParser()
parser.add_argument('--file','-file',help='path to rdb file')
parser.add_argument('--ydata','-ydata',nargs='+',help='column header for the y data (vrad/fwhm/bis_span)')
args = parser.parse_args()
fn = args.file
ydata = args.ydata

# want to ensure this is the rdb file, which contains fwhm columns
assert fn.endswith('.rdb')
# rdb file has rows [header, placeholder, data], so we want to drop the placeholder row
data = pd.read_csv(fn,delimiter='\t')
data = data.drop(index=0)

fig,ax = plt.subplots(len(ydata),sharex='all',squeeze=True,figsize=(5,2*len(ydata)))

for i,var in enumerate(ydata):
    t = np.asarray(data['rjd']).astype(float)
    y = np.asarray(data[var]).astype(float)
    if var=='vrad':
        err = np.asarray(data['svrad']).astype(float)
    elif var=='fwhm':
        err = np.asarray(data['sig_fwhm']).astype(float)
    elif var=='bis_span':
        err = np.asarray(data['sig_bis_span']).astype(float)
    else:
        raise Exception('Unfamiliar ydata header used. Bailing out.')

    f, logp = rv.bgls(t,y,err)

    ax[i].semilogx(1.0/f,logp,label='HARPSN : TOI-1778', linewidth=1.0)
    ax[i].set_xlim(1.0,100.0)
    ax[i].set_ylim(0.0,1.3*np.max(logp))
    ax[i].plot([6.526,6.526],[0.0,100.0],color='black',ls='dashed',lw=0.75)
    ax[i].set_ylabel('log$_{10}$ (Power)')
    ax[i].legend(loc='upper right')

ax[-1].set_xlabel('Period (days)')
plt.tight_layout()
plt.savefig('TOI-1778_BGLS_multivar.png' %ydata,dpi=1000)
plt.show()
