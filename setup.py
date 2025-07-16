import numpy as np
import os,sys
from Feom.hashmap import generateHashmap#,Convert_to_list
from Feom.utils import writeZ,writeI,writeParams,FORT_SWITCHES

npF = np.asfortranarray # Aliasing to make the code more legible

#
#   Setup class for the FEOM integrator
#
class Setup:
    def __init__(self,bath,pot,params):
        ### Add all of the parameters to the class
        self.__dict__.update(vars(params))
        ### Add the bath parameters that are needed
        self.gam_ks = bath.gam_ks
        self.C_ks = bath.C_ks
        self.c_U = bath.c_U
        self.c_D_LEFT = bath.c_D_LEFT
        self.c_D_RIGHT = bath.c_D_RIGHT
        self.lowTcoef = bath.lowTcoef
        self.N_nonmats = bath.N_nonmats
        ### Arguments for the integrator
        self.H_mat = pot.H_mat # Hamiltonian matrix
        self.s_mat = pot.s_mat # perturbation operator (could be q for example) 

        ### Calculate Hashmaps
        # The formatting of the Hashmaps are as follows:
        # I2ind[I] : int I -> list of ints corresponding to the BCF indecies
        # ind2I[ind] :tuple of ints -> int I
        # This is done because lists are not hashable, and tuples are   
        self.ADO_index, self.I0s = generateHashmap(self.K,self.L,self.N_nonmats) #hash map from the index of the ADO to the index of the BCF
        # self.ADO_index, self.I0s = Convert_to_list(self.I2ind) # new indexing for FORTRAN
 

    def generate_input_files(self,x0):
        # Format all of the data
        x0fort =np.zeros((self.Imax,self.ns,self.ns),dtype=complex,order='F')
        for I in range(self.Imax):
            x0fort[I,:,:] = npF(x0[:,:,I])
        # Write the data to the files
        writeZ('Fortrho',x0fort)
        writeI('FortADO_index',npF(self.ADO_index)) 
        writeI('FortI0s',npF(self.I0s)+1) # +1 because fortran is 1 indexed (not 0 indexed like in python)
        writeZ('Fortgam_ks',npF(self.gam_ks))
        writeZ('FortC_ks',npF(self.C_ks))
        writeZ('Fortc_U',npF(self.c_U))
        writeZ('Fortc_D_LEFT',npF(self.c_D_LEFT))
        writeZ('Fortc_D_RIGHT',npF(self.c_D_RIGHT))
        writeZ('FortH_mat',npF(self.H_mat))
        writeZ('Forts_mat',npF(self.s_mat))
        # Write the parameters to the file
        writeParams('Fortparams',self)
        # get the switches and name of the makefile
        makefile_command, self.executable_suffix = FORT_SWITCHES(self)
        self.executable_name=f'propagation{self.executable_suffix}'
        # Copy the correct fortran executable to the temporary directory
        print(f'\n Copying the fortran executable {self.executable_name} to the tmp/ directory\n')
        script_dir = os.path.dirname(os.path.abspath(__file__))
        if os.path.exists(f'{script_dir}/fort/executables/{self.executable_name}'):
            os.system(f' cp {script_dir}/fort/executables/{self.executable_name} ./tmp/') 
        else:
            print(f'\n Need to compile the fortran code. Run the command:\n\n make fast {makefile_command}\n')
            sys.exit()
        return
