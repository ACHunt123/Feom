import numpy as np
import os,sys,glob
import shutil
from pathlib import Path
import Feom.src.fort as fort_pkg
from Feom.src.utils import FORT_SWITCHES,write_sparse,write_Zvec,writeParams,read_Zvec
from Feom.src.initialize.config import SimConfig
import Feom.src.initialize.config_requirements as cfg  
# from Feom.src.hierarchy.real_exps import generate_liouvillian #works for only Debye overdamped
# from Feom.src.hierarchy.complex_exps import generate_liouvillian #works for all BCFs
from Feom.src.hierarchy.z_exps_multibath import generate_liouvillian #works for all BCFs + multiple baths + accelerated

#
#   Setup class for the FEOM integrator with sparse matrices
#
class Setup:
    def __init__(self, config: SimConfig):
        # save the input configuration before we do anything to it
        # config.save('config.json')
        # Unpack the configuration (either inputted through dicts or from a previous json)
        self.bath = self._setup_bath(config.bath)
        self.pot = self._setup_system(config.system)
        self.params = self._setup_params(config.params)
        self._setup_Xi(config.terminator)
        
        ### Generate the sparse matrix Liouvillian
        generate_liouvillian(self)
        print('liouvillian generated')

    def _setup_system(self, sys_dict):
        """
        Derives requirements for the System and builds the object.
        """
        has_s_mat = 's_mat' in sys_dict # single coupling operator (ie single bath)
        has_v_ks = 'V_ks' in sys_dict   # multiple coupling operators (multiple different baths)
        if not (has_s_mat ^ has_v_ks):  # Use XOR operator (^)
            if not has_s_mat and not has_v_ks:
                raise ValueError("System config requires either 's_mat' or 'V_ks'.")
            else:
                raise ValueError("System config cannot have both 's_mat' and 'V_ks'. Pick one!")
        sys_obj= self._create_object_from_dict(
            input_dict=sys_dict,
            required_keys=cfg.REQUIRED_SYS_MANUAL,
            obj_name="ManualPotential")
        # If s_mat specified, all V_ks are the same 
        if has_s_mat:
            sys_obj.V_ks = [sys_obj.s_mat for V_mat in range(len(self.bath.C_ks))]
        # Checks
        self._validate_system(sys_obj)
        # Calculate the derived quantities
        sys_obj.ns=np.shape(sys_obj.H_mat)[0]

        return sys_obj
    
    def _validate_system(self, pot):
        """
        Performs sanity checks on the manually loaded system object.
        """
        for k,V_mat in enumerate(pot.V_ks):
            if np.shape(V_mat) != np.shape(pot.H_mat): # Array lengths
                raise ValueError(f"System Dimension Mismatch: the {k}th coupling operator has shape {np.shape(V_mat)}, but Hamiltonian has shape { np.shape(pot.H_mat)}.")
        try: # Type safety
            pot.V_ks = [np.array(V_mat, dtype=complex) for V_mat in pot.V_ks]
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
        bath_obj.K = len(bath_obj.C_ks)
        # Put in complex conjugate modes in if not
        all_gam_ks = np.unique(np.concatenate([bath_obj.gam_ks, np.conj(bath_obj.gam_ks)]))
        if len(all_gam_ks) != len(bath_obj.gam_ks):
            print(f"adding complex conjugate poles")
            all_C_ks = np.zeros_like(all_gam_ks)
            for k, gam_k in enumerate(all_gam_ks):
                if gam_k in bath_obj.gam_ks:
                    all_C_ks[k] = bath_obj.C_ks[np.where(bath_obj.gam_ks==gam_k)]
            ### Uncomment the below to print the replacements
            # print('OLD')
            # for k, gam_k in enumerate(bath_obj.gam_ks):
            #     print(f'k={k}       {bath_obj.C_ks[k]}-->{bath_obj.gam_ks[k]}\n')
            # print('NEW')
            # for k, gam_k in enumerate(all_gam_ks):
            #     print(f'k={k}       {all_C_ks[k]}-->{all_gam_ks[k]}\n')
            bath_obj.gam_ks = all_gam_ks
            bath_obj.C_ks = all_C_ks
    
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
        params_obj.ns=self.pot.ns
        params_obj.K=self.bath.K

        return params_obj
    
    def _setup_Xi(self, terminator_dict):
        ''' Sets the terminator and gets the flags 
            required given said terminator
        '''   
        if terminator_dict == None: #no terminator added
            self.params.LTCorr = None
            self.Xi = 0+0.j
        else:
            # inherit the attributes and default the normal terminator to be same for each ADO
            self.params.LTCorr = getattr(terminator_dict, 'correction_type', 'same_for_each_ADO') # switch for the compiler flags of FORTRAN

            if  self.params.LTCorr == 'same_for_each_ADO':
                self.Xi = np.zeros((self.params.ns**2,self.params.ns**2),dtype=complex) 
                if np.shape(terminator_dict["Xi"]) != (self.params.ns**2,self.params.ns**2):
                    raise ValueError('Shape of the terminator (one for each ADO) is incorrect')
                self.Xi[:,:]=terminator_dict["Xi"]

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
            works with either the entire set, or the 0th
            vectorizes the rho as well with FORTRAN ordering'''
        x_in_flat=x_in.flatten(order='F')
        self.ADOs = np.zeros(self.params.Ntot, dtype=complex)
        if mode == '0th':  #just set the initial system one
            if len(x_in_flat)!= self.params.ns**2: raise RuntimeError(f"0th ADO has shape {x_in.shape} but should be {self.params.ns}*{self.params.ns}")
            self.ADOs[:len(x_in_flat)] = x_in_flat
        elif mode =='all': #set the entire set 
            if len(x_in_flat)!= self.params.Ntot: raise RuntimeError(f"all ADOs have length {len(x_in_flat)} but should be {self.params.Ntot}")
            self.ADOs = x_in_flat
        else:
            raise ValueError(f"Unknown mode '{mode}'. Use '0th' or 'all'.")
        self._ADOs_loaded=True

    def read_final_ADOs(self):
        ADOs=read_Zvec(f"{self.tmp_folder}/Fortrho.oup", size=(self.params.ns,self.params.ns,self.params.Nados))
        return ADOs


    def generate_input_files(self,tmp_folder='tmp'):
        # Check that everything has been initialized
        if not hasattr(self,"_ADOs_loaded"): raise RuntimeError("ADOs have not been loaded yet")
        
        # Put the files in a temporary folder
        self.tmp_folder = tmp_folder
        self.dest_dir = Path.cwd() / self.tmp_folder
        self.dest_dir.mkdir(parents=True, exist_ok=True) # make the tmp if it doesnt exist

        ## write them all to files
        write_sparse(f"{self.tmp_folder}/FortLiouvillian.inp",self.Liouvillian)
        write_Zvec(f"{self.tmp_folder}/Fortrho.inp", self.ADOs)
        writeParams(f"{self.tmp_folder}/Fortparams.inp", self)

        self._input_files_generated =True


    def insert_executable(self):
        """Fetches the executable installed in the python environment (compiled using pip install .)"""
        if not hasattr(self,"_input_files_generated"): raise RuntimeError("Input files have not been generated yet")
        _, self.executable_suffix = FORT_SWITCHES(self)
        self.executable_name=f'propagation{self.executable_suffix}'
        # Get the directory of the imported fort package
        pkg_dir = Path(fort_pkg.__file__).resolve().parent
        exe_source = pkg_dir / 'executables' / self.executable_name
        # We check the python_environments site-packages directly if the local path fails.
        if not exe_source.exists():
            import site
            # Combine system/venv packages WITH the user packages (~/.local/...)
            all_sites = site.getsitepackages()
            all_sites.append(site.getusersitepackages()) 
            for sp in all_sites:
                candidate = Path(sp) / 'Feom' / 'src' / 'fort' / 'executables' / self.executable_name
                if candidate.exists():
                    exe_source = candidate
                    break       
        if not exe_source.exists():
            print(f"\n[Error] binary '{self.executable_name}' not found.")
            print(f"Verified paths:\n - {pkg_dir}/executables/\n - Environment and User site-packages")
            print("\nTo fix, run: pip install .\n")
            sys.exit(1)
        # copy the executable to the run folder
        dest_file = self.dest_dir / self.executable_name
        shutil.copy2(exe_source, dest_file)
        

    def insert_executable_manually(self):
        """Fetches the executable generated locally via 'make'."""
        if not hasattr(self,"_input_files_generated"): raise RuntimeError("Input files have not been generated yet")
        self.params.LTCorr=None
        makefile_command, self.executable_suffix = FORT_SWITCHES(self)
        self.executable_name=f'propagation{self.executable_suffix}'
        # makefile_command, self.executable_suffix = None,None
        # self.executable_name=f'main'
        # Copy the correct fortran executable to the temporary directory
        print(f'\n Copying the fortran executable {self.executable_name} to the "{self.tmp_folder}" directory\n')
        # get the repo root
        repo_root = Path(__file__).resolve().parent.parent
        # get the executable source code
        exe_source = repo_root / 'fort' / 'executables' / self.executable_name
        # destination directory (tmp folder in current working directory)
        dest_file = self.dest_dir / self.executable_name
        # copy to destination if it exists
        if exe_source.exists():
            shutil.copy2(exe_source, dest_file) 
        else: # Construct the make command 
            print(f'\n[Error] Executable not found at:\n{exe_source}\n')
            print(f'Please compile the fortran code from the repo root:')
            print(f'  {makefile_command}\n')
            sys.exit(1)

    def go(self,extra_commands='',cleanup=True,save_raw=False):
        # Run the executable
        quiet=[' > /dev/null ',''][1] # whether to print out stuff or not
        run_cmd = f'cd {self.tmp_folder}/; {extra_commands} ./{self.executable_name} {quiet}'
        os.system(run_cmd)
        
        #Load the data and format it and attach it to the simulation object
        data= np.loadtxt(f'{self.tmp_folder}/output')
        t= data[:,0]
        ns= self.params.ns
        re_rho = data[:,1:1+ns**2]
        im_rho = data[:,1+ns**2:]
        rho = re_rho + 1.j*im_rho
        self.rho = rho.reshape((len(t),ns,ns), order='C')  # reshape the data to be a 3D array
        self.t_arr = data[:,0]
        raw_labels = ["Time"] + [f"{p}_{i}{j}" for p in ["Re", "Im"] for i in range(ns) for j in range(ns)]
        #Pad every label to be exactly 'w' characters wide (center-aligned)
        formatted_header = "".join([f"{label:^{25}}" for label in raw_labels])
        if save_raw: np.savetxt('raw_output.dat',data,header=formatted_header,fmt='%25.16e', delimiter='')
        #Clean up the temporary directory, and move misc outfiles away
        for file in glob.glob(f'{self.tmp_folder}/*.dat'):
            shutil.move(file, '.')
        if cleanup: self.safe_cleanup() #clean up the temporary directory
        return
    
    def safe_cleanup(self,path=None):
        ''' Safer cleaning up of the tmp directory '''
        if path == None: #default to cleaning up the temporary folder
            path=f"{self.tmp_folder}/"
        if os.path.exists(path) and os.path.isdir(path):
            try:
                shutil.rmtree(path)
                print(f"Successfully cleaned up {path}")
            except OSError as e:
                print(f"Error: {e.strerror}")
        else:
            print(f"Path {path} does not exist or is not a directory.")



    

