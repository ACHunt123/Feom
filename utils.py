import numpy as np
import os


# Some useful functions for the HEOM code



def orthotransmatgen(nb):  #makes the unitary matrix that diagonalises the RPMD hessian
    C = np.zeros((nb,nb),order='F')
    nbo2 = nb//2
    for j in range(0,nb):
        C[j,0] = np.sqrt(1/nb)
        for k in range(1,nbo2):
            C[j,k] = np.sqrt(2/nb)*np.cos((2*j*k*np.pi)/nb)
        C[j,nbo2] = np.sqrt(1/nb)*(-1)**j
        for k in range(nbo2+1,nb):
            C[j,k] = np.sqrt(2/nb)*np.sin((2*j*k*np.pi)/nb)
    return C


def print_progress(current_m1, total, bar_length=40,extra=''):
    current = current_m1+1
    """
    Prints the progress of a task on the same line.

    :param current: Current progress count.
    :param total: Total count for the task.
    :param bar_length: Length of the progress bar.
    """
    progress = current / total
    bar = '#' * int(progress * bar_length) + '-' * (bar_length - int(progress * bar_length))

    string="Progress: ["+bar+"] "+str(current)+"/"+str(total)
    # if current % (total//1) == 0:
    print(string ,end='\r', flush=True)
    if current == total:
        print(' '*len(string),end='\r', flush=True)

### Write functions for the input files in the fortran code
def writeZ(filename, array): # for the large complex arrays (also works for real arrays)
    if not os.path.exists(f"tmp/"): os.makedirs(f"tmp/")
    with open(f"tmp/{filename}", "w") as f:
        f.write("&complex_data\n")
        for idx, c in np.ndenumerate(array): # idx is the index of the array
            ids = ','.join(map(str, [i+1 for i in idx]))
            f.write(f"{c.real:22.15e}\n".replace('e','d'))
            f.write(f"{c.imag:22.15e}\n".replace('e','d'))
            #D22.15 is the format specifier for fortran for each number
        f.write("/\n")
    return
def writeI(filename, array): # for the large integer arrays
    if not os.path.exists(f"tmp/"): os.makedirs(f"tmp/")
    with open(f"tmp/{filename}", "w") as f:
        f.write("&integer_data\n")
        for idx, c in np.ndenumerate(array): # idx is the index of the array
            ids = ','.join(map(str, [i+1 for i in idx]))
            f.write(f"{c:10d}\n")
            #I10 is the format specifier for fortran for each number
        f.write("/\n")
    return
def writeParams(filename,params): # Writes small parameters into file
    if not os.path.exists(f"tmp/"): os.makedirs(f"tmp/")
    with open(f"tmp/{filename}", "w") as f:
        Ktot=params.K+params.N_nonmats
        f.write("Ktot,L,hbar,lowTcoef,Imax,ns,dt,nttot,lowTcoef_switch\n")
        lowTcoef_switch = 1 if hasattr(params, 'LTCorr') else 0
        f.write(f"{Ktot:10d}{params.L:10d}{params.hbar.real:22.15e}{params.lowTcoef.real:22.15e}{params.Imax:10d}{params.ns:10d}{params.dt:22.15e}{params.nttot:10d}{lowTcoef_switch:10d}\n".replace('e','d'))
        f.write("/\n")
    return

# ### Write functions for the output files  and metadata
def out_filename(params): # Name of the output file - contains all the parameters
    namestr = ''
    for key in params.__dict__.keys():
        if key in ['header','out_name','executable_suffix']: continue # dont put the header in the filename
        #get type of variable for formatting
        if type(params.__dict__[key]) == float:
            namestr += f'{key}{params.__dict__[key]:.3g}_'
        if type(params.__dict__[key]) == int:
            namestr += f'{key}{params.__dict__[key]}_'
        if type(params.__dict__[key]) == str:
            namestr += f'{key}{params.__dict__[key]}_'
    ### Clean up the string
    namestr = namestr.replace('bathname','BTH').replace('bathmode','').replace('potname','POT').replace('[N/N]','NoN').replace('[N-1/N]','Nm1oN')
    return f'outfile_{namestr[:-1]}.out'


def FORT_SWITCHES(params):
# Supported compile-time switches (SWITCHES):
#   -DLowTCorr      Enable low-temperature correction via double commutator term
#   -DPrint_ADOs    Print the ADOs to file every N timesteps
#   -DSIA           Use SIA step instead of RK4 step (default)
    switches = []
    if hasattr(params, 'LTCorr') or params.lowTcoef != 0:
        switches.append('LowTCorr')
    if params.print_ADOs:
        switches.append('Print_ADOs')
    if not params.noSIA:
        switches.append('SIA')
    # sort the switches to be alphabetical
    switches.sort()
    makefile_command= f'SWITCHES=" -D{" -D".join(switches)}"'
    executable_suffix='_'+'_'.join(switches)
    if len(switches) == 0:
        executable_suffix = ''
        makefile_command = 'SWITCHES=""'    
    return makefile_command, executable_suffix


def printparams(params):
    ### Make metadata for headers in files - again this contains all the parameters
    runcommand = 'feom.py '
    metadata = "Input parameters used\n"
    metadata += "------------------------------------------------------------------------------------\n"
    for key in params.__dict__.keys():
        if type(params.__dict__[key]) == float:
            metadata += f'{key} = {params.__dict__[key]:.3g} \n'
        if type(params.__dict__[key]) == int:
            metadata += f'{key} = {params.__dict__[key]} \n'
        if type(params.__dict__[key]) == str:
            metadata += f'{key} = {params.__dict__[key]}\n'
        if key not in ['header','out_name', 'executable_suffix']:
            runcommand += f'--{key} {params.__dict__[key]} '
    metadata += 'To compile the FORTRAN executable use the following command (in the Feom/fort directory) \n'
    metadata += f'make fast {FORT_SWITCHES(params)[0]}' + '\n'
    if ('DSIA' in FORT_SWITCHES(params)[0]):
        metadata += 'SIA is Hardcoded with Krylov_dim = 8, Krylov_tol 1e-8, which so far has not had any issues\n'
    metadata += "------------------------------------------------------------------------------------\n"
    metadata += 'To run this code use the following command\n'
    metadata += runcommand + '\n'
    metadata += "------------------------------------------------------------------------------------\n"
    metadata += 'Data \n'
    return metadata
