#!/usr/bin/env python3
from hashmap import  total_length
from utils import print_progress
from integrator import Integrator
import Heom.baths as baths
import Heom.potentials as potentials
from Heom.parser import   params
import numpy as np
import scipy
import sys
import matplotlib.pyplot as plt

if(0):    # Avoid numpy parallelisation
    import os
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
#
# Basic HEOM code
#

### Parameters 
hbar=params.hbar
L = params.L           # the max tier of the ADO expansion
K = params.K           # the number of thermal exponentials in the BCFs - [for debye the number of exponentials = K+1]
ns = params.ns         # number of states to be propagated


### Get bath coefficients
bath = baths.getbath(params.bathname)(params)
C_ks,gam_ks = bath.get_coeffs() 
N_nonmats = bath.N_nonmats # Get the number of non-matsubara terms in the BCF

## Generate matrices in eigenbasis
pot = potentials.getpotential(params.potname)(params)
H_mat = pot.H_matrix()      # Hamiltonian matrix in the eigenbasis
s_mat = pot.pos_matrix()    # position operator matrix in the eigenbasis

### Initial conditions and correlation function function [all hardcoded into potentials]
rho_s0,Zs = pot.initcond() # initial density operator and partition function 
Corr = pot.corr                # function object to calculate the correlation function - also scales time if neccesary

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
Imax = total_length(K,L,N_nonmats)          # the total number of ADOs
rho = np.zeros((ns,ns,Imax),dtype=complex)  # holds all of the ADOs. rho[s,s',0] is the system density matrix
rho[:,:,0] = rho_s0                         # put in the initial density matrix into the ADOS

#setup the integrator
lowTcoef = params.eta/(params.beta*params.hbar**2) - (1/params.hbar**2)*np.sum(np.real(C_ks)/gam_ks) if params.bathmode == 'matsubara' else 0
integrator = Integrator(gam_ks,C_ks,ns,Imax,H_mat,s_mat,K,hbar,L,lowTcoef,N_nonmats)

### Propagation
tmax = 20
dt= 0.001
nt =int(tmax/dt)+1
t_arr = np.arange(nt)*dt


### Plotting
show_plots = 0
if show_plots:
    fig, (ax1, ax2) = plt.subplots(1,2)
    s_arr = np.arange(ns)
    plt.ion()
    if(1):# plot the analytical solution for uncoupled harmonic oscillator
        t_arr,Css_analyt_re = pot.analytic_uncoupled(t_arr=t_arr)
        ax2.plot(t_arr,Css_analyt_re,'--',color='red',label='analytical')
        np.savetxt('CssANALYTIC.txt',np.transpose([t_arr,Css_analyt_re]))
    line, = ax2.plot([],[],'b',label='numerical')


Css = np.zeros_like(t_arr,dtype=complex)
FORTRAN=1
for it in range(nt):

    Css[it], t_arr[it] = Corr(rho[:,:,0],t_arr[it]) #t is an arguement as it may be scaled by potential params
    rho = integrator.rk4_step(rho,dt, FORTRAN=FORTRAN)
    # plot the density matrix and the Css
    if(it%10==0):
        print_progress(it,nt)
        if(show_plots):
            line.remove()
            ax1.clear()
            line, = ax2.plot(t_arr[:it],Css[:it].real,'b',label='numerical')
            ax1.contourf(s_arr,s_arr,np.real(rho[:,:,0]),cmap='viridis')
            plt.pause(0.01)

#write results to file
# NOTE need a file namer
np.savetxt('Css.txt',np.transpose([t_arr,np.real(Css)]))
print('files saved')
plt.plot(t_arr,Css.real)
data=np.loadtxt('CssANALYTIC.txt')
plt.plot(data[:,0],data[:,1],'--')
plt.show()
if show_plots: plt.show()

### Functions


