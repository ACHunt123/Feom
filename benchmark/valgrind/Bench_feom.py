#!/usr/bin/env python3
'''
   +---------------------------------------+
   |   FEOM: Fortran heirarchical          |  
   |       Equations Of Motion             |
   |           By A. C. Hunt 2025          |
   +---------------------------------------+
   THIS FILE IS USED TO BENCHMARK THE DEBYE SPIN-BOSON MODEL
   AGAINST TOM FAY'S MATLAB CODE, AND THUS FORMATS THE OUTPUT 
    TO MATCH HIS OUTPUT.
'''
from Feom.setup import Setup
from Feom.parser import params
import numpy as np


### Load all the parameters into the setup object and generate input files
sim = Setup(params)

### Get the initial conditions from the potential class (ie e^-betaH A from the potential)
rho_s0,Zs = sim.pot.initcond()  # initial density operator and partition function 
rho = np.zeros((params.ns,params.ns,params.Imax),dtype=complex)  # holds all of the ADOs. rho[s,s',0] is the system density matrix
rho[:,:,0] = rho_s0                                              # put in the initial density matrix into the ADOS

### Write the input files
sim.generate_input_files(rho) # writes the inputfiles and fortran code to tmp/ directory


### Go
run_here = 1
if(run_here): # run the fortran code in the temporary directory
    sim.go(extra_commands='valgrind --leak-check=full --show-leak-kinds=all  -s')
