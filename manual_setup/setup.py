import numpy as np
import os,sys
from Feom.hashmap import generateHashmap, total_length
from Feom.baths.utils import get_C_UDs,generate_Terminator
from Feom.utils import writeZ,writeI,writeParams,FORT_SWITCHES,out_filename,printparams

import shutil
from pathlib import Path
npF = np.asfortranarray # Aliasing to make the code more legible

#
#   Setup class for the FEOM integrator
#

import sys
import os
import numpy as np
# from tools import ... 
# import baths

import Feom.manual_setup.config_definitions as cfg  

class ManualSetup:
    def __init__(self, sys_dict, bath_dict, params_dict):
        

        # build the objects from the dictionaries
        self.pot    = self._setup_system(sys_dict)  #calling it pot for now (as the other code uses that)
        self.bath   = self._setup_bath(bath_dict)      
        self.params = self._setup_params(params_dict)  


#         ### Write the parameters to a file and filename
#         params.header = printparams(self)
#         params.out_name = out_filename(self)
        
        ### Calculate the C_U, c_D_LEFT, c_D_RIGHT coefficients for the bath (that are used in the FEOM code)
        get_C_UDs(self.bath,self.params.L)

        ### Calculate Hashmaps for indexing ADOs
        self.ADO_index, self.I0s = generateHashmap(self.bath.K,self.params.L,self.bath.N_nonmats) #hash map from the index of the ADO to the index of the BCF
        self.params.Imax = total_length(self.bath.K,self.params.L,self.bath.N_nonmats)     # the total number of ADOs
        print(f'Total number of ADOs: {self.params.Imax}')

        ### generate the terminator
        # generate_Terminator(self)
        self.Xi=0.0j

    def _setup_system(self, sys_dict):
        """
        Derives requirements for the System and builds the object.
        """
        sys_obj= self._create_object_from_dict(
            input_dict=sys_dict,
            required_keys=cfg.REQUIRED_SYS_MANUAL,
            obj_name="ManualPotential")
        # Checks
        self._validate_system(sys_obj)
        # Calculate the derived quantities
        sys_obj.ns=np.shape(sys_obj.H_mat)[0]

        return sys_obj
    
    def _validate_system(self, pot):
        """
        Performs sanity checks on the manually loaded system object.
        """
        if np.shape(pot.s_mat) != np.shape(pot.H_mat): # Array lengths
            raise ValueError(f"System Dimension Mismatch: coupling operator has shape {np.shape(pot.s_mat)}, but Hamiltonian has shape { np.shape(pot.H_mat)}.")
        try: # Type safety
            pot.s_mat = np.array(pot.s_mat, dtype=complex)
            pot.H_mat = np.array(pot.H_mat, dtype=complex)
        except Exception as e:
            raise TypeError(f"Could not convert system arrays to complex numbers: {e}")

    def _setup_bath(self, bath_dict):
        """
        Derives requirements for the Bath based on its name and builds the object.
        """
        bath_obj = self._create_object_from_dict(input_dict=bath_dict,
            required_keys=cfg.REQUIRED_BATH_MANUAL,
            obj_name="ManualBath")
        # Checks
        self._validate_bath(bath_obj)
        # Calculate derived quantities
        bath_obj.N_exp = len(bath_obj.C_ks)
        bath_obj.N_nonmats = 0
        bath_obj.K = len(bath_obj.C_ks)
    
        return bath_obj
    
    def _validate_bath(self, bath):
        """
        Performs sanity checks on the manually loaded bath object.
        """
        # Array lengths
        if len(bath.C_ks) != len(bath.gam_ks):
            raise ValueError(f"Bath Dimension Mismatch: 'C_ks' has length {len(bath.C_ks)}, but 'gam_ks' has length {len(bath.gam_ks)}.")
        # Type safety
        try:
            bath.C_ks = np.array(bath.C_ks, dtype=complex)
            bath.gam_ks = np.array(bath.gam_ks, dtype=complex)
        except Exception as e:
            raise TypeError(f"Could not convert bath arrays to complex numbers: {e}")
        # Check 3: Zero-Length Check (Edge case)
        if len(bath.C_ks) == 0:
            print("Warning: Bath has 0 modes. Simulation will be purely unitary.")

    def _setup_params(self, params_dict):
        params_obj = self._create_object_from_dict(
            input_dict=params_dict,
            required_keys=cfg.REQUIRED_PARAMS_MANUAL,
            defaults=cfg.DEFAULT_PARAMS_MANUAL,
            obj_name="ManualParams")
        # calculate derived quantities and imports from other objects
        params_obj.hbar=1
        params_obj.ns=self.pot.ns
        params_obj.K=self.bath.K
        # insert compiler defaults

        return params_obj
        
    def _create_object_from_dict(self, input_dict, required_keys, defaults=None, obj_name="GenericObject"):
        """
        Creates an object from a dictionary, applying defaults for missing keys.
        
        Priority:
        1. input_dict (User values)
        2. defaults (Fallback values)
        """
        if defaults is None: defaults = {}
        # check requirements (We only check keys that are NOT in defaults. If a required key has a default, it's safe.)
        effective_keys = set(input_dict.keys()).union(defaults.keys())
        if not required_keys.issubset(effective_keys):
            missing = required_keys - effective_keys
            raise ValueError(f"[{obj_name}] Missing required arguments: {missing}")
        # merge defaults and the data
        final_data = defaults.copy()
        final_data.update(input_dict)
        # create object
        class DynamicContainer:
            def __init__(self, **kwargs):
                for k, v in kwargs.items():
                    setattr(self, k, v)
            def __repr__(self):
                return f"<{obj_name} attributes={list(self.__dict__.keys())}>"

        return DynamicContainer(**final_data)
    
    def generate_input_files(self,x0):
        # Format all of the data
        x0fort =np.zeros((self.params.Imax,self.params.ns,self.params.ns),dtype=complex,order='F')
        for I in range(self.params.Imax):
            x0fort[I,:,:] = npF(x0[:,:,I])
        # Write the data to the files
        writeZ('Fortrho',x0fort)
        writeI('FortADO_index',npF(self.ADO_index)) 
        writeI('FortI0s',npF(self.I0s)+1) # +1 because fortran is 1 indexed (not 0 indexed like in python)
        writeZ('Fortgam_ks',npF(self.bath.gam_ks))
        writeZ('FortC_ks',npF(self.bath.C_ks))
        writeZ('Fortc_U',npF(self.bath.c_U))
        writeZ('Fortc_D_LEFT',npF(self.bath.c_D_LEFT))
        writeZ('Fortc_D_RIGHT',npF(self.bath.c_D_RIGHT))
        writeZ('FortH_mat',npF(self.pot.H_mat))
        writeZ('Forts_mat',npF(self.pot.s_mat))
        writeZ('FortTerminator',npF(self.Xi))
        # Write the parameters to the file
        writeParams('Fortparams',self)
        
    def insert_executable(self):
        # get the switches and name of the makefile

        makefile_command, self.executable_suffix = FORT_SWITCHES(self)
        self.executable_name=f'propagation{self.executable_suffix}'
        # Copy the correct fortran executable to the temporary directory
        print(f'\n Copying the fortran executable {self.executable_name} to the tmp/ directory\n')
        # get the repo root
        repo_root = Path(__file__).resolve().parent.parent  
        # get the executable source code
        exe_source = repo_root / 'fort' / 'executables' / self.executable_name
        # destination directory (tmp folder in current working directory)
        dest_dir = Path.cwd() / 'tmp'
        dest_file = dest_dir / self.executable_name
        # copy to destination if it exists
        if exe_source.exists():
            # shutil is Python's standard tool for copying files
            shutil.copy2(exe_source, dest_file) 
        else: # Construct the make command 
            print(f'\n[Error] Executable not found at:\n{exe_source}\n')
            print(f'Please compile the fortran code from the repo root:')
            print(f'  {makefile_command}\n')
            sys.exit(1)
    

    def go(self,extra_commands='',cleanup=True):
        # Run the executable
        quiet=[' > /dev/null ',''][1]
        os.system(f'cd tmp/; {extra_commands} ./propagation* {quiet}')
        #Load the data and format it
        data= np.loadtxt('tmp/output')
        formatted_data, self.params.header = self.pot.format_output(data,self.params.header)
        #Save it with nice filename and header
        np.savetxt(self.params.out_name,formatted_data.real,header=self.params.header)
        #Clean up the temporary directory
        os.system('mv tmp/*.out .') if os.path.exists('tmp/*.out') else None  # move the output files to the parent directory [only for if we print the ADOs]
        if cleanup: os.system('rm -r tmp/ -f') #clean up the temporary directory
        return 

    

