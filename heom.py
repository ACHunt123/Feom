from hashmap import  total_length
from harm_oscillator import Harmonic_oscillator
from debye_bath import Debye_bath
from integrator import Integrator
import numpy as np
import scipy
import sys
import matplotlib.pyplot as plt


# Avoiding numpy parallelisation
# import os
# os.environ['OPENBLAS_NUM_THREADS'] = '1'
# os.environ['MKL_NUM_THREADS'] = '1'

#
# Basic HEOM code
#


### Parameters
beta600 = 526.2918822622093
beta300 = 1052.5837645244185
beta150 = 2105.167529048837

beta = beta600 #in Adam's code
hbar=1
K = 1           #  the number of elements in the BCFs
L = 0           # the depth of the ADO expansion
ns = 10         # number of eigenstates to be propagated


### Hamiltonian and system setup - harmonic oscillator - same as Adam's code
m =1741.1
d0 = 0.18748
alpha = 1.1605
diff = 2 * d0 *  alpha**2
const = d0 *  alpha**2
omega = (diff/m)**(0.5)

### Position basis parameters for calculation of matrices
xmin = -2
xmax = 15
dx = 0.01


### Bath parameters - Debye bath
bathmode = ['nbead','matsubara'][0]
eta_crit = 2*m*omega  #critical cutoff frequency 
eta=2*eta_crit
gam= omega


### Get bath coefficients
bath = Debye_bath(eta,gam,beta,hbar,K,bathmode)
C_ks,gam_ks = bath.get_coeffs()

## Generate matrices in eigenbasis
pot = Harmonic_oscillator(m,omega,hbar,xmin,xmax,dx,ns)

H_mat = pot.H_matrix()  # Hamiltonian matrix in the eigenbasis
s_mat = pot.pos_matrix() # position operator matrix in the eigenbasis

rho_s0 = pot.rho0(beta)           # initial system density matrix = e^(-beta*H_s) 

#Trace operator
def trace(matrix): return (np.sum(np.diag(matrix)))   

#partition function
Zs = trace(rho_s0)                                    # partition function - calculated without lamb shift



if(0):  # tests
    # Test that the partition function is correct
    theta=hbar*omega*beta
    print(1/(1-np.exp(-theta)),'analytic partition function')
    print(Zs,'partition function calculated','\n')

    # Test that the thermal expectation value of q^2 is correct
    x = np.exp(-beta*hbar*omega)
    print((hbar/(2*m*omega)) * (1+x)/(1-x), 'analytic <q^2>')

    # Calculate the thermal expectation value of q^2 for the system
    print(trace(rho_s0@s_mat@s_mat)/Zs,'<q^2> ') #t
    # When multiplying by the s matrix, we DO NOT ADD in the ds factor because the s matrix would otherwise be divergent
    # lim dx->0[ <x_i|x_j> dx] = delta_ij

### Variables and arrays
Imax = total_length(K,L)                    # the total number of ADOs
rho = np.zeros((ns,ns,Imax),dtype=complex)  # holds all of the ADOs. rho[s,s',0] is the system density matrix
rho[:,:,0] = rho_s0@s_mat                   # initial condition = e^-beta*H_s * s

#setup the integrator
integrator = Integrator(gam_ks,C_ks,ns,Imax,H_mat,s_mat,K,hbar,L)

### Propagation
tmax = 5000
dt= 1
nt =int(tmax/dt)
t_arr = np.linspace(0,tmax,nt)
dt = t_arr[1]-t_arr[0]

### Plotting
fig, (ax1, ax2) = plt.subplots(1,2)
plt.ion()

# plot the analytical solution for uncoupled harmonic oscillator
x = np.exp(beta*hbar*omega/2)
xm1 = x**-1
Css_analyt_re = ((hbar/(2*m*omega)) * (x+xm1)/(x-xm1))*np.cos(omega*t_arr)
ax2.plot(t_arr,Css_analyt_re,'--',color='red',label='analytical')
ax2.set_ylim(-1,1)
Css = np.zeros_like(t_arr)
line, = ax2.plot(t_arr,np.zeros_like(t_arr),'b',label='numerical')

s_arr = np.arange(ns)
for it in range (nt):

    rho = integrator.rk4_step(rho,dt)

    Css[it] = trace(rho[:,:,0]@s_mat)/Zs 

    # plot the density matrix and the Css
    if(it%10==0):
        print(it)
        line.remove()
        ax1.clear()
        line, = ax2.plot(t_arr[:it],Css[:it],'b',label='numerical')
        ax1.contourf(s_arr,s_arr,np.real(rho[:,:,0]),cmap='viridis')
        # print(trace(rho[:,:,0].T@np.diag(s_arr**2)*ds)/np.trace(rho[:,:,0]*ds))
        plt.pause(0.01)
plt.show()

### Functions


