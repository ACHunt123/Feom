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

def writeParams(filename,params,dt,nttot): # Writes small parameters into file
    if not os.path.exists(f"tmp/"): os.makedirs(f"tmp/")
    with open(f"tmp/{filename}", "w") as f:
        f.write("K,L,hbar,lowTcoef,N_nonmats,Imax,ns,dt,nttot\n")
        f.write(f"{params.K:10d}{params.L:10d}{params.hbar.real:22.15e}{params.lowTcoef.real:22.15e}{params.N_nonmats:10d}{params.Imax:10d}{params.ns:10d}{dt:22.15e}{nttot:10d}\n".replace('e','d'))
        f.write("/\n")
    return