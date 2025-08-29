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


    ### Calculate the terminator
    Xi = np.zeros((1,params.ns**2,params.ns**2),dtype=complex)
    I = np.eye(params.ns)
    Vcross = np.kron(pot.s_mat,I) - np.kron(I,pot.s_mat.T)  # commutator superoperator for the system-bath coupling operator
    # Add on the terminator contribution from constant term k in BCF
    k = getattr(bath, 'k', 0)                                           # get the constant term in the BCF if it exists, otherwise 0
    k_term_coef = -2*bath.eta*bath.gam*k/(params.beta*params.hbar**2)   # the constant term in the BCF gives a low temperature correction term
    Xi += k_term_coef * Vcross @ Vcross

    ### Add on the terminator contribution from the terminated frequencies
    LTCorr = getattr(params, 'LTCorr', None)
    # get the Cks and gam_ks for which we have terminated (will not be done explicitly for matsubara as they are infinite)
    bath.gamks_term = bath.gam_ks[bath.N_exp_prop:]
    bath.Cks_term = bath.C_ks[bath.N_exp_prop:]
    if LTCorr is None:
        pass   

    elif LTCorr == 'PT2':   # Second order perturbative terminator on the terminated frequencies
        Ldata = getL(params,pot)
        Xi += get_Xi_n(params,0,Ldata)      # add Xi0 of tom's for each ADO

    elif LTCorr == 'IT':    # Markovian (Ishizaki-Tanimura) terminator on the terminated frequencies
        if bath.mode == 'matsubara':        # The IT Term. coeff. has infinite mats terms in it
            IT_term_coef = bath.eta* ((1/(2*bath.hbar))* ((1/(np.tan(bath.beta*bath.hbar*bath.gam/2))) - (2/(bath.beta*bath.hbar*bath.gam))))  ### Terms without removing of the matsubara terms that have been included
            IT_term_coef -= bath.eta*(2*bath.gam/(bath.beta*bath.hbar**2))*np.sum(1/(bath.gam**2*np.ones_like(bath.ws[1:]) - bath.ws[1:]**2))  ### remove the Matsubara terms that have been explicitly included
            Xi += IT_term_coef * Vcross @ Vcross
        else:
            Xi += get_IT(params,pot,bath)   # Calculate the IT terminator from the unused C_ks and gam_ks (calculated in Bath class)

    elif LTCorr == 'NZ2':   # Tom Fay's Nakajima-Zwanwig terminator on the terminated frequencies
        Xi_n = np.zeros((params.Imax,params.ns**2,params.ns**2),dtype=complex)  # one xi for EACH ADO (will overwrite Xi)
        Ldata = getL(params,pot)
        for I in range(params.Imax):
            Xi_n[I,:,:] += get_Xi_n(params,I,Ldata)
            Xi_n[I,:,:] += Xi               # add on the contributions from Xi onto each Xi_n
        Xi = Xi_n   # override Xi to be the full set of ADO-specific terminators

    else:
        raise ValueError('Invalid LTCorr type')
    
    ### Assign a dummy variable if the terminator is zero
    if (abs(Xi)==0).all():
        print('No terminator added')
        sim.Xi = np.zeros((1,1,1),dtype=complex) # this is just a dummy var
    else:
        sim.Xi = Xi
    return
        

    # Calculate the C_ks and gam_ks for a given set of ws
    def get_Xi_n(self,n):
        ''' Get the Xi_n matrix for a given ADO n (list of indices)'''

        return 





