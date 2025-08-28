import numpy as np
from types import SimpleNamespace
import sys

def get_C_UDs(params):
    ''' Generate the coefficients C_U, c_D_LEFT, c_D_RIGHT from the C_ks for the bath (that are used in the FEOM code)'''
    # Calculate the coefficients C_U, c_D_LEFT, c_D_RIGHT for the bath (that are used in the FEOM code)
    params.c_U = np.zeros((params.N_exp,params.L+1),dtype=complex)
    params.c_D_LEFT = np.zeros((params.N_exp,params.L+1),dtype=complex)
    params.c_D_RIGHT = np.zeros((params.N_exp,params.L+1),dtype=complex)
    for ki in range(params.N_exp):
        for nk in range(params.L+1):
            params.c_U[ki,nk] = np.sqrt((nk+1)*abs(params.C_ks[ki]))
            if abs(params.C_ks[ki]) < 1e-10:
                params.c_D_LEFT[ki,nk] = 0.0
                params.c_D_RIGHT[ki,nk] = 0.0
                print(f'Warning: C_ks({ki}) is zero, setting superoperator terms to zero')
            else:
                params.c_D_LEFT[ki,nk] = -np.sqrt(nk/abs(params.C_ks[ki]))*params.C_ks[ki]
                params.c_D_RIGHT[ki,nk] = np.sqrt(nk/abs(params.C_ks[ki]))*np.conj(params.C_ks[ki])
    return 


def generate_Terminator(sim):
    '''Function to generate terminator for HEOM.
    IT - ishizaki-Tanimura terminator (same for each ADO)
    PT2 - 2nd order perturbative terminator (same for each ADO)
    NZ2 - Nakajima-Zwanwig terminator (different for each ADO, more expensive)


    self.init_lowtcoef is the low-temp correction either due to IT, or the k term from AAA and Pade[N/N] it will be added to
    It is a class, as it need to store its type and size
    '''
    params = sim.params
    pot= sim.pot
    bath = sim.bath

    ### Calculate the terminator contribution from low temperature correction
    I = np.eye(params.ns)
    Vcross = np.kron(pot.s_mat,I) - np.kron(I,pot.s_mat.T)  # commutator superoperator for the system-bath coupling operator
    Xi_lowtcorr =  bath.lowTcoef * Vcross @ Vcross  # low temperature correction term (same for all ADOs)

    ### Get the Liouvillian and transformation matrices (for the NZ2 terminator)
    def getL(params,pot):
        ''' Get the Liouvillian matrix for the system Hamiltonian, and transformation matrices'''
        I = np.eye(params.ns)
        Lsys = -1.j*(np.kron(pot.H_mat,I) - np.kron(I,pot.H_mat.T))/params.hbar
        eigvals, eigvecs = np.linalg.eig(Lsys)
        Lams = np.diag(eigvals)
        Pis = eigvecs
        Pis_inv = np.linalg.inv(Pis)
        if(0):#test the eigen-decomposition
            L_reconstructed = Pis @ Lams @ Pis_inv
            error = np.linalg.norm(Lsys - L_reconstructed)
            print(f"Eigen-decomposition reconstruction error: {error:.2e}")
        return SimpleNamespace(Lsys=Lsys,Lams=np.diag(eigvals),Pis=eigvecs,Pis_inv=np.linalg.inv(eigvecs))



    def get_Xi_n(params,n,Ldata):
        ''' Get the Xi_n matrix for a given ADO n (list of indices)'''
        sys.exit('Not implemented yet')

    def get_IT(params,pot,bath):
        ''' Get the IT terminator matrix, calculated from the unused C_ks and gam_ks in the bath object'''
        I = np.eye(params.ns)
        VL = np.kron(pot.s_mat,I)
        VR = np.kron(I,pot.s_mat.T)
        Vcross = VL-VR  # commutator superoperator for the system-bath coupling operator
        Xi_IT = np.zeros((params.ns**2,params.ns**2),dtype=complex)
        nterm=len(bath.gamks_term)
        for k in range(nterm):
            gamk = bath.gamks_term[k]
            Ck = bath.Cks_term[k]
            O=(Ck*VL -Ck.conj()*VR)
            Xi_IT -= Vcross@O/(gamk*params.hbar**2)
        return Xi_IT



    ### Add on the terminator contribution from the terminated frequencies
    LTCorr = getattr(params, 'LTCorr', None)
    #get the Cks and gam_ks for which we have terminated
    bath.gamks_term = bath.gam_ks[bath.N_exp_prop:]
    bath.Cks_term = bath.C_ks[bath.N_exp_prop:]

    if LTCorr is None:      # no extra low temperature correction
        if bath.lowTcoef == 0: 
            sim.Xi = np.zeros((1,1,1),dtype=complex) # this is just a dummy var
        else:
            sim.Xi = np.zeros((1,params.ns**2,params.ns**2),dtype=complex) 
            sim.Xi[0,:,:] = Xi_lowtcorr  # add on the low temperature correction term from constant k
        return

    elif LTCorr == 'PT2': # second order perturbative terminator on the terminated frequencies
        sim.Xi = np.zeros((1,params.ns,params.ns),dtype=complex)
        Ldata = getL(params,pot)
        params.Xi[0,:,:] = get_Xi_n(params,0,Ldata) # add Xi0 of tom's for each ADO
        params.Xi[0,:,:] += Xi_lowtcorr             # add on the low temperature correction term from constant k

    elif LTCorr == 'IT': # Markovian (Ishizaki-Tanimura-like) terminator on the terminated frequencies
        sim.Xi = np.zeros((1,params.ns**2,params.ns**2),dtype=complex)
        sim.Xi[0,:,:] = get_IT(params,pot,bath)      # add the markovian terms the terminated frequencies
        sim.Xi[0,:,:] += Xi_lowtcorr        # add on the low temperature correction term from constant k


    elif LTCorr == 'NZ2': # Tom Fay's Nakajima-Zwanwig terminator on the terminated frequencies
        sim.Xi = np.zeros((params.Imax,params.ns**2,params.ns**2),dtype=complex) # one xi for EACH ADO
        Ldata = getL(params,pot)
        for I in range(params.Imax):
            sim.Xi[I,:,:] = get_Xi_n(params,I,Ldata)
            sim.Xi[I,:,:] += Xi_lowtcorr    # add on the low temperature correction term from constant k

    else:
        raise ValueError('Invalid LTCorr type')
    

    return
    sys.exit()

        

    # Calculate the C_ks and gam_ks for a given set of ws
    def get_Xi_n(self,n):
        ''' Get the Xi_n matrix for a given ADO n (list of indices)'''

        return 





