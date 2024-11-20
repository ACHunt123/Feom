from hashmap import  total_length
from harm_oscillator import H_matrix, rho0
from integrator import Integrator
import numpy as np
import scipy
import matplotlib.pyplot as plt
#
# Basic HEOM code
#


### Parameters
beta = 0.1
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



# rho_s0 = scipy.linalg.expm(-beta*H_mat) #THIS DOESN'T WORK NOT SURE WHY

H_mat, x_arr, dx = H_matrix(ns,m,omega,hbar,nx,xmin,xmax)  # Hamiltonian matrix
s_mat = np.diag(x_arr)                           # the system operator
rho_s0 = rho0(beta,ns,m,omega,hbar,nx,xmin,xmax)           # initial system density matrix

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

# Iprop = #numbers of ADOs to propagete on the current time step



rho = np.zeros((ns,ns,Imax),dtype=complex)  # holds all of the ADOs. rho[s,s',0] is the system density matrix
rho[:,:,0] = rho_s0                         # initial condition

integrator = Integrator(ns,Imax,H_mat,s_mat,K,hbar,max_N)

### Propagation
dt=0.01
fig = plt.figure()
ax = fig.add_subplot(111)
plt.ion()
for i in range (0,1000):

    rho = integrator.rk4_step(rho,dt)

    # plot the density matrix
    if(i%10==0):

        ax.contourf(x_arr,x_arr,np.real(rho[:,:,0]),cmap='viridis')
        print(np.trace(rho[:,:,0].T@np.diag(x_arr**2)*dx)/np.trace(rho[:,:,0]*dx))
        plt.pause(0.01)
plt.show()

### Functions


