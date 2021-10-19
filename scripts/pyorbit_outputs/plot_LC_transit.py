import math
import numpy as np
import matplotlib.pyplot as plt

LC_data_full = np.genfromtxt('LCdata_full.dat',skip_header=1)
LC_data_lc_model_b = np.genfromtxt('LCdata_lc_model_b.dat',skip_header=1)

#plt.plot(LC_data_full[:,0],LC_data_full[:,1])
plt.scatter(LC_data_lc_model_b[:,2],LC_data_lc_model_b[:,3],s=8,color='black')
plt.scatter(LC_data_lc_model_b[:,2],LC_data_lc_model_b[:,7],s=6,color='red')

plt.ylim(0.998,1.0015)
plt.xlim(0.19,0.24)

plt.savefig('LC_transit.png')
plt.show()
