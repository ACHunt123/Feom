#!/usr/bin/env python3
# File: heom.py
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

beta = beta150 #in Adam's code
hbar=1
L = 3           # the depth of the ADO expansion
K = 3           #  the number of elements in the BCFs
ns = 10         # number of eigenstates to be propagated


### Hamiltonian and system setup - harmonic oscillator - same as Adam's code
m =1741.1
d0 = 0.18748
alpha = 1.1605
diff = 2 * d0 *  alpha**2
const = d0 *  alpha**2
omega = (diff/m)**(0.5)


### Position basis parameters for calculation of matrices
xmin = -5
xmax = 5
dx = 0.01


### Bath parameters - Debye bath
bathmode = ['nbead','matsubara'][1]
eta_crit = 2*m*omega  #critical cutoff frequency 
eta_ADAM=2*eta_crit
eta = eta_ADAM*omega 
gam= omega


### Get bath coefficients
bath = Debye_bath(eta,gam,beta,hbar,K,bathmode)
C_ks,gam_ks = bath.get_coeffs() 

## Generate matrices in eigenbasis
pot = Harmonic_oscillator(m,omega,hbar,xmin,xmax,dx,ns)

H_mat = pot.H_matrix()  # Hamiltonian matrix in the eigenbasis
s_mat = pot.pos_matrix() # position operator matrix in the eigenbasis
rho_s = pot.rho0(beta)  #unrormed system density matrix = e^(-beta*H_s) 

Zs = np.trace(rho_s)     #partition function
rho_s0 = rho_s@s_mat/Zs  # initial density operator, normalised by the partition function 


if(0):  # tests
    print('### Tests ###')
    # Test that the partition function is correct
    theta=hbar*omega*beta
    print(1/(1-np.exp(-theta)),'analytic partition function')
    print(Zs,'partition function calculated','\n')

    # Test that the thermal expectation value of q^2 is correct
    x = np.exp(-beta*hbar*omega)
    print((hbar/(2*m*omega)) * (1+x)/(1-x), 'analytic <q^2>')

    # Calculate the thermal expectation value of q^2 for the system
    print(np.trace(rho_s0@s_mat),'<q^2> ') #t
    # When multiplying by the s matrix, we DO NOT ADD in the ds factor because the s matrix would otherwise be divergent
    # lim dx->0[ <x_i|x_j> dx] = delta_ij
    print('')

### Variables and arrays
Imax = total_length(K,L)                    # the total number of ADOs
rho = np.zeros((ns,ns,Imax),dtype=complex)  # holds all of the ADOs. rho[s,s',0] is the system density matrix
rho[:,:,0] = rho_s0                         # put in the initial density matrix into the ADOS

#setup the integrator
lowTcoef = eta/(beta*hbar**2) - (1/hbar**2)*np.sum(np.real(C_ks)/gam_ks) if bathmode == 'matsubara' else 0
integrator = Integrator(gam_ks,C_ks,ns,Imax,H_mat,s_mat,K,hbar,L,lowTcoef)

### Propagation
tmax = 200
dt= 1
nt =int(tmax/dt)+1
t_arr = np.arange(nt)*dt


### Plotting
show_plots = 0
if show_plots:
    fig, (ax1, ax2) = plt.subplots(1,2)
    plt.ion()


if(0):# plot the analytical solution for uncoupled harmonic oscillator
    x = np.exp(beta*hbar*omega/2)
    xm1 = x**-1
    Css_analyt_re = ((hbar/(2*m*omega)) * (x+xm1)/(x-xm1))*np.cos(omega*t_arr)
    ax2.plot(t_arr,Css_analyt_re,'--',color='red',label='analytical') if show_plots else [None]
    # ax2.set_ylim(-0.02,0.02)
    line, = ax2.plot(t_arr,np.zeros_like(t_arr),'b',label='numerical') if show_plots else [None]
    np.savetxt('CssANALYTIC.txt',np.transpose([t_arr,Css_analyt_re]))


Css = np.zeros_like(t_arr,dtype=complex)
s_arr = np.arange(ns)

for it in range(nt):

    Css[it] = np.trace(rho[:,:,0]@s_mat)
    rho = integrator.rk4_step(rho,dt)


    # plot the density matrix and the Css
    if(it%10==0 and show_plots):
        line.remove()
        ax1.clear()
        line, = ax2.plot(t_arr[:it],Css[:it],'b',label='numerical')
        ax1.contourf(s_arr,s_arr,np.real(rho[:,:,0]),cmap='viridis')
        plt.pause(0.01)

#write results to file
np.savetxt('Css.txt',np.transpose([t_arr,np.real(Css)]))
print('files saved')
print(Css.shape)
plt.show()

### Functions


