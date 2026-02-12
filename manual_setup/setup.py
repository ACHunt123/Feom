import numpy as np
import os,sys,glob
from Feom.utils.hashmap import generateHashmap, total_length
from Feom.baths.utils import get_C_UDs,generate_Terminator
from Feom.utils.utils import writeZ,writeI,writeParams,FORT_SWITCHES,out_filename,printparams

import shutil
from pathlib import Path
npF = np.asfortranarray # Aliasing to make the code more legible

#
#   Setup class for the FEOM integrator
#

import sys
import os
import numpy as np
from Feom.manual_setup.config import SimConfig
# from tools import ... 
# import baths

import Feom.manual_setup.config_requirements as cfg  

class ManualSetup:
    def __init__(self, config: SimConfig):
        # save the input configuration before we do anything to it
        config.save('config.json')
        # Unpack the configuration (either inputted through dicts or from a previous json)
        self.pot = self._setup_system(config.system)
        self.bath = self._setup_bath(config.bath)
        self.params = self._setup_params(config.params)
        
        # Calculate the C_U, c_D_LEFT, c_D_RIGHT coefficients for the bath (that are used in the FEOM code)
        get_C_UDs(self.bath,self.params.L)

        # Calculate Hashmaps for indexing ADOs
        self.ADO_index, self.I0s = generateHashmap(self.bath.K,self.params.L,self.bath.N_nonmats) #hash map from the index of the ADO to the index of the BCF
        self.params.Imax = total_length(self.bath.K,self.params.L,self.bath.N_nonmats)     # the total number of ADOs
        print(f'Total number of ADOs: {self.params.Imax}')
        
        ### generate the terminator
        # generate_Terminator(self)
        self._setup_Xi(config.terminator)

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
    
    def _setup_Xi(self, terminator_dict):
        ''' Sets the terminator and gets the flags 
            required given said terminator
        '''   
        if terminator_dict == None: #no terminator added
            self.params.LTCorr = None
            self.Xi = np.zeros((1,1,1),dtype=complex)
        else:
            # inherit the attributes and default the normal terminator to be same for each ADO
            self.params.LTCorr = getattr(terminator_dict, 'correction_type', 'same_for_each_ADO') # switch for the compiler flags of FORTRAN

            if  self.params.LTCorr == 'same_for_each_ADO':
                self.Xi = np.zeros((1,self.params.ns**2,self.params.ns**2),dtype=complex) 
                if np.shape(terminator_dict["Xi"]) != (self.params.ns**2,self.params.ns**2):
                    raise ValueError('Shape of the terminator (one for each ADO) is incorrect')
                self.Xi[0,:,:]=terminator_dict["Xi"]

            elif  self.params.LTCorr == 'different_for_each_ADO':
                self.Xi = np.zeros((self.params.Imax,self.params.ns**2,self.params.ns**2),dtype=complex) 
                if np.shape(terminator_dict["Xi"]) != (self.params.Imax,self.params.ns**2,self.params.ns**2):
                    raise ValueError('Shape of the terminator (different for each ADO) is incorrect')
                self.Xi[0,:,:]=terminator_dict["Xi"]
            
            else:
                raise ValueError('correction_type should be in [same_for_each_ADO,different_for_each_ADO]') 

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
    
    def set_initial_ADOs(self,x_in,mode='0th'):
        ''' Setup the initial state of the ADOs
            works with either the entire set, or the 0th'''
        self.x0 = np.zeros((self.params.ns, self.params.ns, self.params.Imax), dtype=complex)
        if mode == '0th': #just set the initial system one
            self.x0[:,:,0]=x_in
        elif mode =='all': #set the entire set 
            self.x0=x_in
        else:
            raise ValueError(f"Unknown mode '{mode}'. Use '0th' or 'all'.")
        self._ADOs_loaded=True

    def generate_input_files(self):
        # Check that everything has been initialized
        if not hasattr(self,"_ADOs_loaded"): raise RuntimeError("ADOs have not been loaded yet")
        # Format all of the data
        x0fort =np.zeros((self.params.Imax,self.params.ns,self.params.ns),dtype=complex,order='F')
        for I in range(self.params.Imax):
            x0fort[I,:,:] = npF(self.x0[:,:,I])
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
            shutil.copy2(exe_source, dest_file) 
        else: # Construct the make command 
            print(f'\n[Error] Executable not found at:\n{exe_source}\n')
            print(f'Please compile the fortran code from the repo root:')
            print(f'  {makefile_command}\n')
            sys.exit(1)
    
    def go(self,extra_commands='',cleanup=True,save_raw=True):
        # Run the executable
        quiet=[' > /dev/null ',''][1]
        os.system(f'cd tmp/; {extra_commands} ./propagation* {quiet}')
        #Load the data and format it and attach it to the simulation object
        data= np.loadtxt('tmp/output')
        t= data[:,0]
        ns= self.params.ns
        re_rho = data[:,1:1+ns**2]
        im_rho = data[:,1+ns**2:]
        rho = re_rho + 1.j*im_rho
        self.rho = rho.reshape((len(t),ns,ns), order='C')  # reshape the data to be a 3D array
        self.t_arr = data[:,0]
        header = "Time " + " ".join([f"{part}_{i}{j}" for part in ['Re', 'Im'] for j in range(ns) for i in range(ns)])
        raw_labels = ["Time"] + [f"{p}_{i}{j}" for p in ["Re", "Im"] for i in range(ns) for j in range(ns)]
        #Pad every label to be exactly 'w' characters wide (center-aligned)
        formatted_header = "".join([f"{label:^{25}}" for label in raw_labels])
        if save_raw: np.savetxt('raw_output.dat',data,header=formatted_header,fmt='%25.16e', delimiter='')
        #Clean up the temporary directory, and move misc outfiles away
        for file in glob.glob('tmp/*.dat'):
            shutil.move(file, '.')
        if cleanup: self._safe_cleanup("tmp/") #clean up the temporary directory
        return 
    
    def _safe_cleanup(self,path):
        ''' Safer cleaning up of the tmp directory '''
        if os.path.exists(path) and os.path.isdir(path):
            try:
                shutil.rmtree(path)
                print(f"Successfully cleaned up {path}")
            except OSError as e:
                print(f"Error: {e.strerror}")
        else:
            print(f"Path {path} does not exist or is not a directory.")



    

