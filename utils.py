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
        f.write("Ktot,L,hbar,lowTcoef,Imax,ns,dt,nttot\n")
        f.write(f"{Ktot:10d}{params.L:10d}{params.hbar.real:22.15e}{params.lowTcoef.real:22.15e}{params.Imax:10d}{params.ns:10d}{params.dt:22.15e}{params.nttot:10d}\n".replace('e','d'))
        f.write("/\n")
    return

# ### Write functions for the output files  and metadata
# def out_filename(potkey,params): #no potkey added for the moment
#     if mode in ['MQCL','QMdensitymat','QMwavefunc_splitop']:
#         return f'{mode}_Nx{Nx}_dt{dt*fs:.3g}_tmax{tmax*fs:.3g}_m{m:.1g}_xa{xa:2g}_xb{xb:.2g}.txt'.replace(' ','') 
#     elif mode in ['QMwavefunc_DVR']:
#         return f'{mode}_Nx{Nx}_Ne{Ne}_dt{dt*fs:.3g}_tmax{tmax*fs:.3g}_m{m:.1g}_xa{xa:2g}_xb{xb:.2g}.txt'.replace(' ','')
#     else:
#         print('Invalid mode')
#         sys.exit()

# def TM_printparams(potkey,potname,mode,Nx,tmax,m,xa,xb,dt=None,Ne=None):
#     fname = 'parameters' # make parameters file
#     # if not int(input('Do you want to save the parameters used in a file? (1/0)')): return
#     with open(fname,'w') as f:
#         f.write('Input parameters used\n')
#         f.write('---------------------\n')
#         # all parameters
#         f.write(f'potkey = {potkey}\n')
#         f.write(f'mode = {mode}\n')
#         f.write(f'Nx = {Nx}\n')
#         f.write(f'tmax = {tmax*fs}fs\n')
#         f.write(f'm = {m}\n')
#         f.write(f'xa = {xa}\n')
#         f.write(f'xb = {xb}\n')
#         f.write(f'dt = {dt*fs}fs\n\n')
#         # parameters specific to modules
#         if mode == 'QMwavefunc_DVR':
#             f.write(f'Ne = {Ne}\n')
#         # the argparse to generate said parameters
#         f.write('To generate the same parameters, use the following command in target directory (use -batch for multiple to ensure plotting doesnt happen):\n')
#         if mode in ['QMwavefunc_DVR']:
#             f.write(f'wfTLS_DVR.py -Nx {Nx} -Ne {Ne} -xa {xa} -xb {xb} -dt {dt*fs:4g} -tmax {tmax*fs:4g} -m {m} -potname {potname}\n')
#         if mode in ['QMwavefunc_splitop']:
#             f.write(f'wfTLS_splitop.py -Nx {Nx} -xa {xa} -xb {xb} -dt {dt*fs:4g} -tmax {tmax*fs:4g} -m {m} -potname {potname}\n')
#         if mode in ['MQCL','QMdensitymat']:
#             f.write(f'TruncMoyalTLS.py -Nx {Nx} -xa {xa} -xb {xb} -dt {dt*fs:4g} -tmax {tmax*fs:4g} -m {m} -simulation {mode} -potname {potname}\n')
#         f.write('Note: The rest of the parameters are hardcoded for now, so you will need to change them manually\n')
#     ### Make metadata for headers in files
#     metadata = "Input parameters used\n"
#     metadata += "------------------------------------------------------------------------------------\n"
#     metadata += f'potkey = {potkey}\n'
#     metadata += f'mode = {mode}\n'
#     metadata += f'Nx = {Nx}\n'
#     metadata += f'tmax = {tmax * fs} fs\n'
#     metadata += f'm = {m}\n'
#     metadata += f'xa = {xa}\n'
#     metadata += f'xb = {xb}\n'
#     if dt is not None:
#         metadata += f'dt = {dt * fs} fs\n\n'
#     # Parameters specific to the module
#     if mode == 'QMwavefunc_DVR' and Ne is not None:
#         metadata += f'Ne = {Ne}\n'
#     # Command to reproduce the parameters in the target directory
#     if mode != None: metadata += '\n To generate the same parameters, use the following command in the target directory (use -batch for multiple to ensure plotting does not happen):\n'
#     # Generate the appropriate command based on the mode
#     if mode == 'QMwavefunc_DVR':
#         metadata += f'wfTLS_DVR.py -Nx {Nx} -Ne {Ne} -xa {xa} -xb {xb} -dt {dt * fs:4g} -tmax {tmax * fs:4g} -m {m} -potname {potname}\n'
#     elif mode == 'QMwavefunc_splitop':
#         metadata += f'wfTLS_splitop.py -Nx {Nx} -xa {xa} -xb {xb} -dt {dt * fs:4g} -tmax {tmax * fs:4g} -m {m} -potname {potname}\n'
#     elif mode in ['MQCL', 'QMdensitymat']:
#         metadata += f'TruncMoyalTLS.py -Nx {Nx} -xa {xa} -xb {xb} -dt {dt * fs:4g} -tmax {tmax * fs:4g} -m {m} -simulation {mode} -potname {potname}\n'
#     # Add additional note for hardcoded parameters
#     if mode != None: metadata += '\n Note: The rest of the parameters are hardcoded for now, so you will need to change them manually\n'
#     metadata += "------------------------------------------------------------------------------------\n"
#     if mode != None: metadata += 'Data'

#     return metadata
