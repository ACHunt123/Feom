import numpy as np
import matplotlib.pyplot as plt
import sys

'''Function to generate terminator for HEOM.
IT - ishizaki-Tanimura terminator (same for each ADO)
PT2 - 2nd order perturbative terminator (same for each ADO)
NZ2 - Nakajima-Zwanwig terminator (different for each ADO, more expensive)


self.init_lowtcoef is the low-temp correction either due to IT, or the k term from AAA and Pade[N/N] it will be added to
It is a class, as it need to store its type and size
'''
class Terminator():
    def __init__(self,params):
        # Bathmode and settings
        self.init_lowtcoef = params.lowTcoef 
        self.H=params.H_mat
        self.s=params.s_mat
        self.ns=params.ns
        self.hbar=params.hbar
        self.getL()
        sys.exit()

    def getL(self):
        ''' Get the Liouvillian matrix for the system Hamiltonian, and transformation matrices'''
        I = np.eye(self.ns)
        self.Lsys = -1.j*(np.kron(self.H,I) - np.kron(I,self.H.T))/self.hbar
        eigvals, eigvecs = np.linalg.eig(self.Lsys)
        self.Lams = np.diag(eigvals)
        self.Pis = eigvecs
        self.Pis_inv = np.linalg.inv(self.Pis)
        if(0):#test the eigen-decomposition
            L_reconstructed = self.Pis @ self.Lams @ self.Pis_inv
            error = np.linalg.norm(self.Lsys - L_reconstructed)
            print(f"Eigen-decomposition reconstruction error: {error:.2e}")
        return 
        

    # Calculate the C_ks and gam_ks for a given set of ws
    def get_Xi_n(self,n):
        ''' Get the Xi_n matrix for a given ADO n (list of indices)'''

        return 

    # output TCF for a given set of C_ks and gam_ks
    def TCF(self,plotme=False,ax=plt,mode=None):
        if mode is None: mode = self.mode #allowing override of the mode from the __init__

        self.get_coeffs(mode=mode)

        t = np.linspace(0,5,1000)
        C = np.zeros_like(t,dtype=complex)
        # C += self.C_ks[0]*np.exp(-self.gam_ks[0]*t)
        print(f'C0: {self.C_ks[0]}, mode : {mode}')
        for k in range(1,self.N_exp):
            C += self.C_ks[k]*np.exp(-self.gam_ks[k]*t)
        if plotme:
            ax.plot(t,C.real)
            print(mode)
            plt.show()
        return t,C

    # Calculate the C_ks and gam_ks for a bath, mode is the bath decomposition mode
    def get_coeffs(self,mode=None):
        if mode is None: mode = self.mode #allowing override of the mode from the __init__

        if mode == 'highT':
            return np.array([self.C0hot]),np.array([self.gam]) #the high temperature limit - no matsubara terms

        if mode in ['matsubara','nmats']: # Generate the K matsubara frequencies (the nmats will have the same freqs, but different c0 and no low temp truncation)
            betaN = self.beta/self.N_mds
            wN=1/(betaN*self.hbar)
            wns = np.array([2*wN*np.pi*k/self.N_mds for k in range(0,self.mu+1)])
            print(f'mode {mode} with {self.mu} pairs of matsubara modes')
            # print(wns)
            return self.calc_coefs(wns)

        if mode == 'nbead':
            betaN =  self.beta/self.N_mds
            wN=1/(betaN*self.hbar)
            wks = np.array([2*wN*np.sin(np.pi*k/self.N_mds) for k in range(0,self.mu+1)])
            print(f'mode {mode} with {self.mu} pairs of matsubara modes')
            # print(wks)
            return self.calc_coefs(wks)
        else:
            raise ValueError('Invalid mode')
        return

    def get_C_UDs(self):
        # Calculate the coefficients C_U, c_D_LEFT, c_D_RIGHT for the bath (that are used in the FEOM code)
        self.c_U = np.zeros((self.N_exp,self.L+1),dtype=complex)
        self.c_D_LEFT = np.zeros((self.N_exp,self.L+1),dtype=complex)
        self.c_D_RIGHT = np.zeros((self.N_exp,self.L+1),dtype=complex)
        for ki in range(self.N_exp):
            for nk in range(self.L+1):
                self.c_U[ki,nk] = np.sqrt((nk+1)*abs(self.C_ks[ki]))
                if abs(self.C_ks[ki]) < 1e-10:
                    self.c_D_LEFT[ki,nk] = 0.0
                    self.c_D_RIGHT[ki,nk] = 0.0
                    print(f'Warning: C_ks({ki}) is zero, setting superoperator terms to zero')
                else:
                    self.c_D_LEFT[ki,nk] = -np.sqrt(nk/abs(self.C_ks[ki]))*self.C_ks[ki]
                    self.c_D_RIGHT[ki,nk] = np.sqrt(nk/abs(self.C_ks[ki]))*np.conj(self.C_ks[ki])
        return 



