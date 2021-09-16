import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-file')
parser.add_argument('-nplanets')
args = parser.parse_args()
file = args.file
npl = int(args.nplanets)

results = np.genfromtxt(file,dtype=str,delimiter=',sjh124')

# select relevant section of results file to take values from
lstart = [i for i,line in enumerate(results) if line=='Statistics on the physical parameters obtained from the posteriors samples'][0]
lend = [i for i,line in enumerate(results) if 'Parameters corresponding to the Maximum a Posteriori probability' in line][0]

cropped_results = results[lstart:lend]

P = [line for line in cropped_results if line.startswith('P ')]
a = [line for line in cropped_results if line.startswith('a_AU_(M) ')]
e = [line for line in cropped_results if line.startswith('e ')]
i = [line for line in cropped_results if line.startswith('i ')]
K = [line for line in cropped_results if line.startswith('K ')]
M = [line for line in cropped_results if line.startswith('M_Me ')]
R = [line for line in cropped_results if line.startswith('R_Re ')]
Tc = [line for line in cropped_results if line.startswith('Tc ')]
Prot = [line for line in cropped_results if line.startswith('Prot ')]
Pdec = [line for line in cropped_results if line.startswith('Pdec ')]
Oamp = [line for line in cropped_results if line.startswith('Oamp ')]

planet_ids = ['b','c','d','e','f','g']
planet_dict = {}
for j in range(npl):
    planet_dict[planet_ids[j]] = {
        'P': P[j].split()[1:4],
        'a': a[j].split()[1:4],
        'e': e[j].split()[1:4],
        'i': i[j].split()[1:4],
        'K': K[j].split()[1:4],
        'M': M[j].split()[1:4],
        'R': R[j].split()[1:4],
        'Tc': Tc[j].split()[1:4]
    }
print(planet_dict)

if Prot:
    activity_dict = {
        'Prot' : Prot[0].split()[1:4],
        'Pdec' : Pdec[0].split()[1:4],
        'Oamp' : Oamp[0].split()[1:4]
    }
    print(activity_dict)
