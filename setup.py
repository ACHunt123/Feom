import numpy as np
import os,sys
from Feom.hashmap import generateHashmap, total_length
from Feom.baths.utils import get_C_UDs,generate_Terminator
from Feom.utils import writeZ,writeI,writeParams,FORT_SWITCHES,out_filename,printparams
import Feom.baths as baths
import Feom.potentials as potentials

npF = np.asfortranarray # Aliasing to make the code more legible

#
#   Setup class for the FEOM integrator
#
class Setup:
    def __init__(self,params):
        ### Add all of the generalparameters to the class
        self.params=params

        ### Setup bath and potential 
        self.bath = baths.getbath(self.params.bathname)(self.params)
        self.pot = potentials.getpotential(self.params.potname)(self.params)

        ### Write the parameters to a file and filename
        params.header = printparams(self)
        params.out_name = out_filename(self)
        
        ### Calculate the C_U, c_D_LEFT, c_D_RIGHT coefficients for the bath (that are used in the FEOM code)
        get_C_UDs(self.bath)

        ### Calculate Hashmaps
        self.ADO_index, self.I0s = generateHashmap(self.params.K,self.params.L,self.bath.N_nonmats) #hash map from the index of the ADO to the index of the BCF
        self.params.Imax = total_length(self.params.K,self.params.L,self.bath.N_nonmats)     # the total number of ADOs

        ### Calculate the terminator if needed, and add to self
        generate_Terminator(self)


    def generate_input_files(self,x0):
        # Format all of the data
        x0fort =np.zeros((self.params.Imax,self.params.ns,self.params.ns),dtype=complex,order='F')
        for I in range(self.params.Imax):
            x0fort[I,:,:] = npF(x0[:,:,I])
        # Write the data to the files
        writeZ('Fortrho',x0fort)
        writeI('FortADO_index',npF(self.ADO_index)) 
        writeI('FortI0s',npF(self.I0s)+1) # +1 because fortran is 1 indexed (not 0 indexed like in python)
        writeZ('Fortgam_ks',npF(self.bath.gam_ks[:self.bath.N_exp_prop]))
        writeZ('FortC_ks',npF(self.bath.C_ks[:self.bath.N_exp_prop]))
        writeZ('Fortc_U',npF(self.bath.c_U))
        writeZ('Fortc_D_LEFT',npF(self.bath.c_D_LEFT))
        writeZ('Fortc_D_RIGHT',npF(self.bath.c_D_RIGHT))
        writeZ('FortH_mat',npF(self.pot.H_mat))
        writeZ('Forts_mat',npF(self.pot.s_mat))
        writeZ('FortTerminator',npF(self.Xi))
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
    
    def go(self,extra_commands=''):
        # Run the executable
        os.system(f'cd tmp/; {extra_commands} ./propagation*')
        #Load the data and format it
        data= np.loadtxt('tmp/output')
        formatted_data, self.params.header = self.pot.format_output(data,self.params.header)
        #Save it with nice filename and header
        np.savetxt(self.params.out_name,formatted_data.real,header=self.params.header)
        #Clean up the temporary directory
        os.system('mv tmp/*.out .') if os.path.exists('tmp/*.out') else None  # move the output files to the parent directory [only for if we print the ADOs]
        os.system('rm -r tmp/ -f') #clean up the temporary directory
        return 
