import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-nplanet','--nplanet')
parser.add_argument('-ndata','--ndata')
args = parser.parse_args()
npl = int(args.nplanet) # number of planets modelled
ndata = int(args.ndata) # number of data files

# planet labels for plot titles
plabels = ['planet_e','planet_d','planet_c','planet_b'][::-1]
# data labels for plot legends
dlabels = ['HiRES','HARPS-N'][::-1]

for i in range(ndata):
    if ndata == 1:
        planetdat = ['RV_radial_velocities_e.dat','RV_radial_velocities_d.dat',
                     'RV_radial_velocities_c.dat','RV_radial_velocities_b.dat'][::-1][:npl]
        modeldat = ['RV_planet_e_pha.dat','RV_planet_d_pha.dat','RV_planet_c_pha.dat',
                    'RV_planet_b_pha.dat'][::-1][:npl]
    elif ndata > 1:
        planetdat = ['RV%i_radial_velocities_e.dat' %(i+1), 'RV%i_radial_velocities_d.dat' %(i+1),
                     'RV%i_radial_velocities_c.dat' %(i+1), 'RV%i_radial_velocities_b.dat' %(i+1)][::-1][:npl]
        modeldat = ['RV_planet_e_pha.dat','RV_planet_d_pha.dat','RV_planet_c_pha.dat',
                    'RV_planet_b_pha.dat'][::-1][:npl]
    for j, dfile in enumerate(planetdat):
        data = np.genfromtxt(dfile, skip_header = 1)

        print(np.shape(data)[0])

        plt.figure((i*npl)+j,figsize=(6,5))
        gs = gridspec.GridSpec(2,1,height_ratios=[3,1])

        plt.subplot(gs[0])
        plt.scatter(data[:,2],data[:,8],color='black',s=15,marker='o',label=dlabels[i])
        plt.errorbar(data[:,2],data[:,8],yerr=data[:,9],linestyle='None',color='black',elinewidth=0.5)

        for k in range(0,np.shape(data)[0]):
          print(k)
          if (data[k,2]-1.0) > -0.5 and data[k,2]-1.0 < 0.0:
            plt.scatter(data[k,2]-1.0,data[k,8],s=15,marker='o',color='grey')
            plt.errorbar(data[k,2]-1.0,data[k,8],yerr=data[k,9],linestyle='None',color='grey',elinewidth=0.5)
          if data[k,2] + 1.0 > 1.0 and data[k,2]+1.0 < 1.5:
            plt.scatter(data[k,2]+1.0,data[k,8],s=15,marker='o',color='grey')
            plt.errorbar(data[k,2]+1.0,data[k,8],yerr=data[k,9],linestyle='None',color='grey',elinewidth=0.5)

        data_model = np.genfromtxt(modeldat[j],skip_header=1)

        plt.plot(data_model[:,0],data_model[:,1])

        plt.xlim(-0.5,1.5)
        plt.ylim(-15.0,15.0)

        plt.legend(loc='upper right',fontsize=8)

        plt.ylabel('RV [ms$^{-1}$]')
        plt.title(plabels[j])
        #plt.text(-0.45, 25.0, '$K_b$ = 5.90 $\pm$ 1.13 m s$^{-1}$', fontsize = 10)
        #plt.text(-0.45, 20.00, '$M_b$ = 58.1 $\pm$ 11.5 M$_\oplus$', fontsize = 10)
        #plt.text(-0.45, 10.0, '$e_b$ = 0 (assumed)', fontsize = 8)
        #plt.text(1.4,-13.5,'(a)',fontsize=10)

        plt.subplot(gs[1])
        plt.scatter(data[:,2],data[:,10],s=15,marker='o',color='black')
        plt.errorbar(data[:,2],data[:,10],yerr=data[:,11],linestyle='None',color='black',elinewidth=0.5)

        for k in range(0,np.shape(data)[0]):
          if (data[k,2]-1.0) > -0.5 and data[k,2]-1.0 < 0.0:
            plt.scatter(data[k,2]-1.0,data[k,10],s=15,marker='o',color='grey')
            plt.errorbar(data[k,2]-1.0,data[k,10],yerr=data[k,11],linestyle='None',color='grey',elinewidth=0.5)
          if data[k,2] + 1.0 > 1.0 and data[k,2]+1.0 < 1.5:
            plt.scatter(data[k,2]+1.0,data[k,10],s=15,marker='o',color='grey')
            plt.errorbar(data[k,2]+1.0,data[k,10],yerr=data[k,11],linestyle='None',color='grey',elinewidth=0.5)

        plt.plot([-0.5,1.5],[0.0,0.0],color='black')

        plt.xlim(-0.5,1.5)
        plt.ylim(-15.0,15.0)

        plt.xlabel('Phase')
        plt.ylabel('O-C [ms$^{-1}$]')

        plt.savefig('RV_phase_%s_%s.png' %(plabels[j], dlabels[i]))
plt.show()
