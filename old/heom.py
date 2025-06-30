#!/usr/bin/env python3
'''
   +---------------------------------------+
   |   FEOM: Fortran heirarchical          |  
   |       Equations Of Motion             |
   |           By A. C. Hunt 2025          |
   +---------------------------------------+

   This implementation can calculate correlation functions
   with 3 different methods: 
    FULLFORTRAN : 1 to generate input files and run the fortran code
    FORTRAN     : 1 to propagate with f2py vvstep
    else        : propagate with python
'''
from hashmap import  total_length
from utils import print_progress
from integrator import Integrator
import Feom.baths as baths
import Feom.potentials as potentials
from Feom.parser import   params
import numpy as np
import scipy
import sys,os
import matplotlib.pyplot as plt

if(0):    # Avoid numpy parallelisation
    import os
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
#
# Basic HEOM code
#
### Switches
FORTRAN=1       # 1 to propagate with f2py vvstep, 0 to propagate with python
FULLFORTRAN=1   # 1 do generate input files and run the fortran code

### Get bath coefficients
bath = baths.getbath(params.bathname)(params)

## Generate matrices in eigenbasis
pot = potentials.getpotential(params.potname)(params)

### Initial conditions and correlation function function [all hardcoded into potentials]
rho_s0,Zs = pot.initcond()  # initial density operator and partition function 
Corr = pot.corr             # function object to calculate the correlation function - also scales time if neccesary


### Variables and arrays
params.Imax = total_length(params.K,params.L,bath.N_nonmats)          # the total number of ADOs
rho = np.zeros((params.ns,params.ns,params.Imax),dtype=complex)  # holds all of the ADOs. rho[s,s',0] is the system density matrix
rho[:,:,0] = rho_s0                         # put in the initial density matrix into the ADOS


### Integrator
integrator = Integrator(bath,pot,params)

### Plotting
show_plots = 0
if show_plots:
    fig, (ax1, ax2) = plt.subplots(1,2)
    s_arr = np.arange(params.ns)
    plt.ion()
    if(1):# plot the analytical solution for uncoupled harmonic oscillator
        t_arr,Css_analyt_re = pot.analytic_uncoupled(t_arr=params.t_arr)
        ax2.plot(t_arr,Css_analyt_re,'--',color='red',label='analytical')
        np.savetxt('CssANALYTIC.txt',np.transpose([t_arr,Css_analyt_re]))
    line, = ax2.plot([],[],'b',label='numerical')


### Run the fortran code (if neccesary)
if(FULLFORTRAN):
    print('### Generating input files for fortran ###')
    integrator.generate_input_files(rho)
    run_here = True
    if(run_here): # run the fortran code in the temporary directory
        print('### Running fortran code ###')
        os.system('cd tmp/; ./propagation')
        # Load the data from the fortran execution and save it to a file with header showing params
        fname = header = 'Css.txt'#TM_out_filename(potkey,simulation,Nx,dt,tmax,m,xa,xb)
        data= np.loadtxt('tmp/output')
        np.savetxt(fname,data,header=header)
        # os.system(f'mv tmp/output {fname}')
        # os.system('rm -r tmp/ -f') #clean up the temporary directory
    else:
        print('Files and executables ready to go' )
    sys.exit()

### Propagate the system    
Css = np.zeros_like(params.t_arr,dtype=complex)
for it in range(params.nttot):

    Css[it], params.t_arr[it] = Corr(rho[:,:,0],params.t_arr[it]) #t is an arguement as it may be scaled by potential params
    rho = integrator.rk4_step(rho, FORTRAN=FORTRAN) # Propagate either with python or fortran
    # plot the density matrix and the Css
    if(it%10==0):
        print_progress(it,params.nttot)
        if(show_plots):
            line.remove()
            ax1.clear()
            line, = ax2.plot(params.t_arr[:it],Css[:it].real,'b',label='numerical')
            ax1.contourf(s_arr,s_arr,np.real(rho[:,:,0]),cmap='viridis')
            plt.pause(0.01)

#write results to file
# NOTE need a file namer
np.savetxt('Css.txt',np.transpose([params.t_arr,np.real(Css)]))
print('files saved')
plt.plot(params.t_arr,Css.real)
data=np.loadtxt('CssANALYTIC.txt')
plt.plot(data[:,0],data[:,1],'--')
plt.show()
if show_plots: plt.show()

### Functions


