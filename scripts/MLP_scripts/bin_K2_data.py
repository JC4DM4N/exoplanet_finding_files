import numpy as np
import math
import matplotlib.pyplot as plt

outfile = open('Kepler-103_notransits_binned.dat','w')

data = np.genfromtxt('output13.dat')

data_binarray = [[],[],[],[],[],[]]

time_array = []
K2_array = []

nbins = np.shape(data)[0]

for i in range(0,nbins):
  time_binned = 0.0
  data_binned = 0.0
  stdev_binned = 0.0
  stdev = 0.0

  k = 0
  for j in range(0,10):
    if (i*10 + j <= nbins-9):
      time_binned = time_binned + data[i*10 + j,0]
      data_binned = data_binned + data[i*10 + j,1]
      stdev_binned = stdev_binned + data[i*10 + j,2]
      stdev = stdev + ((data[i*10+j,1]-np.mean(data[i*10:i*10+10,1]))*(data[i*10+j,1]-np.mean(data[i*10:i*10+10,1])))
      k = k + 1

  if (time_binned > 0.0):
    time_binned = time_binned/k
    data_binned = data_binned/k
    stdev_binned = stdev_binned/k
    stdev = math.sqrt(stdev/k)

    time_array.append(time_binned)
    K2_array.append(data_binned)

    print(time_binned, data_binned, stdev_binned, 2, 2, -1)

  data_stdev = np.std(data[i*5:i*5+5,1])

  data_binarray.append([time_binned,data_binned,data_stdev,0,0,0])

plt.plot(time_array, K2_array)
plt.plot(data[:,0],data[:,1],ls='dotted')

plt.show()
