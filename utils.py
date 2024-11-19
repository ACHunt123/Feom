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