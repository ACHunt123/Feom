import numpy as np



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