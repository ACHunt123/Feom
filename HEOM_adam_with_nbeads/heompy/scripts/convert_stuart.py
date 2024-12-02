#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from matplotlib import animation

print('Plotting script starts.')

file_dir = os.getcwd() + '/' # + '/output/'

file_name0 = sys.argv[1]
file_path0 = file_dir + file_name0
datafile = np.genfromtxt(file_path0)#, skip_header = 2)
N_s0 = 51
N_t0 = int(datafile.shape[0]/N_s0)


ss0 = np.zeros(N_s0)
ts0 = np.zeros(N_t0)
data0_raw = np.zeros([N_t0+1,N_s0+1])
for s in range(N_s0):
    ss0[s] = datafile[s,1]
    data0_raw[0, s+1] = ss0[s]

for t in range(N_t0):
    ts0[t] = datafile[N_s0*t,0]
    #data[t+1,0] = ts[t]
    data0_raw[t+1,0] = ts0[t]
    for s in range(N_s0):
        data0_raw[t+1,s+1] = datafile[t*N_s0+s,2]
for item in ss0:
    print(item)

for item in ts0:
    print(item)

with open("bench.dat", "w") as outf:
    #for line in data0_raw:
    np.savetxt(outf, data0_raw, fmt='%30.20E')
    #np.savetxt(outf, data0_raw, fmt='%f')

for z in range(len(data0_raw[:])):
    continue

outf.close()
