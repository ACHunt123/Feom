#!/usr/bin/env python3
'''
   +---------------------------------------+
   |   FEOM: Fortran heirarchical          |  
   |       Equations Of Motion             |
   |           By A. C. Hunt 2025          |
   +---------------------------------------+
'''
from hashmap import  total_length
from utils import print_progress,printparams
from setup import Setup
import Feom.baths as baths
import Feom.potentials as potentials
from Feom.parser import   params
import numpy as np
import scipy
import sys,os
import matplotlib.pyplot as plt


### Setup bath, potential and integrator
bath = baths.getbath(params.bathname)(params)
pot = potentials.getpotential(params.potname)(params)

### Initial conditions and correlation function function [all hardcoded into potentials]
rho_s0,Zs = pot.initcond()  # initial density operator and partition function 

### Setup initial density matrix (direct product of system and bath)
params.Imax = total_length(params.K,params.L,bath.N_nonmats)     # the total number of ADOs
rho = np.zeros((params.ns,params.ns,params.Imax),dtype=complex)  # holds all of the ADOs. rho[s,s',0] is the system density matrix
rho[:,:,0] = rho_s0                                              # put in the initial density matrix into the ADOS

### Load all the parameters into the setup object and generate input files
setup = Setup(bath,pot,params)
setup.generate_input_files(rho) # writes the inputfiles and fortran code to tmp/ directory

### Go
run_here = 1
if(run_here): # run the fortran code in the temporary directory
   os.system(f'cd tmp/; ./{setup.executable_name}')
   fname = header = 'Css.txt'
   # get the raw output file and then format it, updating the header to get data_labels
   data= np.loadtxt('tmp/output')
   formatted_data, params.header = pot.format_output(data,params.header)

   np.savetxt(setup.out_name,formatted_data.real,header=params.header)
   os.system('mv tmp/*.out .') if os.path.exists('tmp/*.out') else None  # move the output files to the parent directory [only for if we print the ADOs]
   os.system('rm -r tmp/ -f') #clean up the temporary directory


