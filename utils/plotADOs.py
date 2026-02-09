#!/usr/bin/env python3
'''   A code to plot the ADOs from the FEOM code
    '''

# Firstly, we need to import the ADO indices from ADO_index.out

import numpy as np
import matplotlib.pyplot as plt
import os,sys

indices = np.loadtxt('ADO_index.out',dtype=int) # this is formatted as I then k_1, k_2 ...k_K

# Now we need to format the list such that there
Imax= indices.shape[0]
tiers = np.zeros((Imax),dtype=int) # this will hold the tier and the number of ADOs in that tier 
for I in range(Imax):
    k = indices[I,1:] # this is the k vector for the I-th ADO
    tiers[I] = np.sum(k) # the tier is the sum of the k vector

print('Tiers:',tiers)
# Now open the ados
ados = np.loadtxt('ADOs.out') # this is formatted as I then the ADOs
t = ados[:,0] # the first column is the time
ados = ados[:,1:] # the rest are the ADOs

# setup an interactive plot for the ADOs
# plt.ion()  # Turn on interactive mode
plt.figure(figsize=(10, 6))
for I in range(Imax):
    plt.plot(t, ados[:, I], label=f'ADO {I} (Tier {indices[I,1:]})')
plt.xlabel('Time')
plt.ylabel('ADO Value')
plt.title('ADOs over Time')
plt.legend()
plt.grid()
plt.show()
# Print the shape of the ADOs
print('ADOS shape:',ados.shape)