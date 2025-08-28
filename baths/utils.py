import numpy as np
from types import SimpleNamespace


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
    def getL(params):
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
    
    def get_Xi_n(params,n):
        ''' Get the Xi_n matrix for a given ADO n (list of indices)'''


    LTCorr = getattr(params, 'LTCorr', None) # what type of terminator are we using?
    Ldata = getL(params)                     # get the Liouvillian data   
    I = np.eye(params.ns)

    Vcross = np.kron(pot.s_mat,I) - np.kron(I,pot.s_mat.T)  # commutator superoperator for the system-bath coupling operator
    Xi_lowtcorr=  bath.lowTcoef * Vcross @ Vcross  # low temperature correction term (same for all ADOs)


    sim.Xi = np.zeros((1,params.ns**2,params.ns**2),dtype=complex)
    sim.Xi[0,:,:] = Xi_lowtcorr

    # if LTCorr is None and params.lowTcoef == 0: # no terminator and no constant term in 
    #     params.Xi = np.zeros((1,1,1),dtype=complex) 
    #     return
    
    # elif LTCorr == 'NZ2':
    #     params.Xi = np.zeros((params.Imax,params.ns,params.ns),dtype=complex) # one xi for each ADO
    #     for n in range(params.Imax):
    #         params.Xi[n,:,:] = get_Xi_n(params,n,Ldata)


    # elif LTCorr in ['IT','PT2'] or params.lowTcoef!=0:
    #     params.Xi = np.zeros((1,params.ns,params.ns),dtype=complex) # only one xi for all ADOs
    #     params.Xi[n,:,:] = get_Xi_n(params,0,Ldata)
    #     # fist add on the lowtcoef term

    #     #then add on the Xi_0 term to all of them

    # else:
    #     raise ValueError('Invalid LTCorr type')
    

    return
    sys.exit()

        

    # Calculate the C_ks and gam_ks for a given set of ws
    def get_Xi_n(self,n):
        ''' Get the Xi_n matrix for a given ADO n (list of indices)'''

        return 





