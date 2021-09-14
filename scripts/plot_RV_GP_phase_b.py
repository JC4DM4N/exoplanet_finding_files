import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

data_HARPSN = np.genfromtxt('RV_radial_velocities_b.dat', skip_header = 1)

print(np.shape(data_HARPSN)[0])

plt.figure(figsize=(8,6))
gs = gridspec.GridSpec(2,1,height_ratios=[3,1])

#plt.subplot(gs[0])

#plt.scatter(data_HARPSN[:,0],data_HARPSN[:,3]-data_HARPSN[:,5],color='black',s=15,marker='o',label='HARPS-N')
#plt.errorbar(data_HARPSN[:,0],data_HARPSN[:,3]-data_HARPSN[:,5],yerr=data_HARPSN[:,4],linestyle='None',color='black',elinewidth=0.5)

data_model_kep = np.genfromtxt('RV_radial_velocities_b_full.dat',skip_header=1)

#plt.plot(data_model_kep[:,0],data_model_kep[:,3])

data_model_gp = np.genfromtxt('RV_gp_quasiperiodic_full.dat',skip_header=1)

plt.plot(data_model_gp[:,0],data_model_gp[:,2]+data_model_kep[:,2],zorder=1)

plt.fill_between(data_model_gp[:,0],(data_model_gp[:,2]+data_model_kep[:,2])-data_model_gp[:,3],(data_model_gp[:,2]+data_model_kep[:,2])+data_model_gp[:,3],alpha=0.5)

plt.scatter(data_HARPSN[:,0],data_HARPSN[:,3]-data_HARPSN[:,5],color='black',s=15,marker='o',label='HARPS-N',zorder=2)
plt.errorbar(data_HARPSN[:,0],data_HARPSN[:,3]-data_HARPSN[:,5],yerr=data_HARPSN[:,4],linestyle='None',color='black',elinewidth=0.5,zorder=2)

plt.xlim(9175.0,9360.0)
#plt.ylim(-15.0,15.0)

#plt.legend(loc='upper right',fontsize=8)

plt.ylabel('RV [ms$^{-1}$]')
#plt.title('Planet b')
#plt.text(-0.45, 13.0, '$K_b$ = 4.66 $\pm$ 0.76 m s$^{-1}$', fontsize = 9)
#plt.text(-0.45, 11.00, '$M_b$ = 10.05 $\pm$ 1.46 M$_\oplus$', fontsize = 9)
#plt.text(-0.45, 9.0, '$e_b$ = 0.74 $\pm$ 0.05', fontsize = 9)

#plt.text(1.4,-13.5,'(a)',fontsize=10)

#plt.subplot(gs[1])
#plt.scatter(data_HARPSN[:,2],data_HARPSN[:,10],s=15,marker='o',color='black')
#plt.errorbar(data_HARPSN[:,2],data_HARPSN[:,10],yerr=data_HARPSN[:,11],linestyle='None',color='black',elinewidth=0.5)

#for i in range(0,np.shape(data_HARPSN)[0]):
#  if (data_HARPSN[i,2]-1.0) > -0.5 and data_HARPSN[i,2]-1.0 < 0.0:
#    plt.scatter(data_HARPSN[i,2]-1.0,data_HARPSN[i,10],s=15,marker='o',color='grey')
#    plt.errorbar(data_HARPSN[i,2]-1.0,data_HARPSN[i,10],yerr=data_HARPSN[i,11],linestyle='None',color='grey',elinewidth=0.5)
#  if data_HARPSN[i,2] + 1.0 > 1.0 and data_HARPSN[i,2]+1.0 < 1.5:
#    plt.scatter(data_HARPSN[i,2]+1.0,data_HARPSN[i,10],s=15,marker='o',color='grey')
#    plt.errorbar(data_HARPSN[i,2]+1.0,data_HARPSN[i,10],yerr=data_HARPSN[i,11],linestyle='None',color='grey',elinewidth=0.5)

#plt.plot([-0.5,1.5],[0.0,0.0],color='black')

#plt.xlim(-0.5,1.5)
#plt.ylim(-15.0,15.0)

plt.xlabel('Date (BJD - 2450000)')
#plt.ylabel('O-C [ms$^{-1}$]')

#plt.savefig('RV_phase_b_GP_TOI1778.png')
plt.show()
