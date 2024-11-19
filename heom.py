from hashmap import generateHashmap, total_length
from harm_oscillator import H_matrix, rho0
import numpy as np
import scipy
import matplotlib.pyplot as plt
#
# Basic HEOM code
#


### Parameters
beta = 1
ns = 50         # number of states/xtics (same thing really as we are in the position basis)
K = 3           #  the number of elements in the BCFs
max_N = 2       # the maximum number of elements in the ADOs

### Bath parameters - Debye bath
gam_ks = np.ones(K,dtype=complex)           # the gammas for the BCFs
C_ks = np.ones(K,dtype=complex)             # the prefactors for the BCFs

### Hamiltonian and system setup [tbc]
m =1
omega = 1
hbar = 1

nx = ns
xmin = -10
xmax = 10

# Hamiltonian matrix and initial conditions
H_mat,x_arr,dx = H_matrix(ns,m,omega,hbar,nx,xmin,xmax)

# rho_s0 = scipy.linalg.expm(-beta*H_mat) #THIS DOESN'T WORK NOT SURE WHY
rho_s0 = rho0(beta,ns,m,omega,hbar,nx,xmin,xmax)  # initial system density matrix

if(0):  # tests
    # Test that the thermal expectation value of q^2 is correct
    x = np.exp(-beta*hbar*omega)
    print((hbar/(2*m*omega)) * (1+x)/(1-x))
    print(np.trace(rho_s0.T@np.diag(x_arr**2)*dx)/np.trace(rho_s0*dx))

    # Test that the thermal density matrix is stationary
    print(((rho_s0@H_mat - H_mat@rho_s0)*dx).max(),'maximal deviation from stationarity')




### Variables and arrays
# I             : the natural numbers that are used to index the ADOs
# I2ind         : hash map from the index of the ADO to the index of the BCF in tuple form [e.g. '0,0,0' for the 0th BCF with K=3]
# ind2I         : the inverse map of I2ind
Imax = total_length(K,max_N)                # the total number of ADOs
I2ind , ind2I = generateHashmap(K,max_N)
# The formatting of the Hashmaps are as follows:
# I2ind[I] : int I -> list of ints corresponding to the BCF indecies
# ind2I[ind] :tuple of ints -> int I
# This is done because lists are not hashable, and tuples are


rho = np.zeros((ns,ns,Imax),dtype=complex)  # holds all of the ADOs. rho[s,s',0] is the system density matrix
rho[:,:,0] = rho_s0                         # initial condition

### Functions

# returns the index of the ADO that is one element different from the input ADO
def I_nk_plusminus(I,k,pm):
    ind = I2ind[I]
    nkpm = ind[k]+pm
    # if nkpm is negative, the ADO does not exist and the function returns -1
    if nkpm < 0:
        return -1
    # Otherwise, the function find the index of the ADO with the same elements as ind, but with the kth element changed to nkpm and return it
    else:
        ind[k] = nkpm
        ind = str(ind).replace('[','').replace(']','').replace(' ','')
        return ind2I[ind]  

# gives the unperturbed Liouvillian acting on the ADO = i/hbar [H,rho]
def L0(rho):
    return 1.j/hbar * (H_mat@rho - rho@H_mat)


def drhodt():

    gradient = np.zeros((ns,ns,Imax),dtype=complex)

    for I in range(0,Imax):
        n_ks = np.array(I2ind[I]) # the indexes of the ADO - will be used as coefficients 

        gradient[:,:,I] = -L0(rho[:,:,I]) + np.sum(n_ks*gam_ks)*rho[:,:,I]

        for k in range(K):
            # the ADOs that are one element different from the current ADO
            I_nkp1 = I_nk_plusminus(I,k,+1)
            I_nkm1 = I_nk_plusminus(I,k,-1)

            # The gradient for the +1 terms
            gradient[:,:,I] += 1.j/hbar * commutator(S,rho[:,:,I_nkp1]) * np.sqrt(C_ks[k]*(n_ks[k]+1)) ###redo#

            # The gradient for the -1 terms [may not exist]
            if I_nkm1 != -1:
                gradient[:,:,I] += 1.j/hbar * commutator(S,rho[:,:,I_nkp1]) * np.sqrt(C_ks[k]*(n_ks[k]+1)) ###redo#

    return gradient

