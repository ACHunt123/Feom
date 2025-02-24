import numpy as np
import sys  
from scipy import special
import matplotlib.pyplot as plt


### Generates the Hamiltonian matrix spin boson model in diabatic representation
# Note adiabatic repn is without _'s and diabatic repn is with _'s on matrices
class Spin_boson:
    def __init__(self,ns=2):
        assert ns == 2, "Spin boson model is hardcoded for two levels"
        self.ns = ns
        self.Delta = 1
        self.eps = 1*self.Delta
        # compute the diabatic<->adiabatic transformation matrix ONCE
        H_mat = np.zeros((self.ns,self.ns),dtype=complex)
        H_mat[0,0] = -self.eps
        H_mat[1,1] = self.eps
        H_mat[1,0] = self.Delta
        H_mat[0,1] = self.Delta
        self.H_mat = H_mat
        self.eigs,self.Uda = np.linalg.eigh(H_mat) # diagonalize the Hamiltonian

    # Generates the Hamiltonian matrix in the diabatic basis
    def H_matrix(self):
        return self.H_mat

    def pos_matrix(self): # perturbation matrix in diabatic basis
        s_mat = np.zeros((self.ns,self.ns),dtype=complex)
        s_mat[0,0] = -1
        s_mat[1,1] = 1
        return s_mat

    # Generates the initial system density matrix and the partition function
    def initcond(self,beta): #(here we just put the system in the excited state)
        rho_s0 = np.zeros((self.ns,self.ns),dtype=complex)
        rho_s0[1,1] = 1
        Zs = 1
        return rho_s0,Zs

    # Generate the correlation function (this is hardcoded)
    def corr(self,rho_s,t): # (return population of the excited diabatic state)
        return rho_s[1,1],t

    def analytic_uncoupled(self,beta,t_arr=np.arange(0,10,0.1)):# calculatees the analytic solution for the uncoupled system
        phi_0=np.zeros((self.ns),dtype=complex)
        phi_0[1]=1
        phi_0 = self.Uda.T@phi_0 #transform into eigenbasis
        pop = np.zeros_like(t_arr)
        for i,t in enumerate(t_arr):
            K = np.diag(np.exp(-1j*self.eigs*t))
            phi_t = K@phi_0
            pop[i] = np.abs((self.Uda@phi_t)[1])**2
        return t_arr,pop