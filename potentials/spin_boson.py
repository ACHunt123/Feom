import numpy as np
import sys  
from scipy import special
import matplotlib.pyplot as plt


### Generates the Hamiltonian matrix spin boson model in diabatic representation
class Spin_boson:
    def __init__(self,ns=2):
        assert ns == 2, "Spin boson model is hardcoded for two levels"
        self.ns = ns
        self.eps = 1
        self.Delta = 1
        self.m = 1

    # Generates the Hamiltonian matrix in its eigenbasis
    def H_matrix(self):
        Hmat = np.zeros((self.ns,self.ns),dtype=complex)
        Hmat[0,0] = -self.eps
        Hmat[1,1] = self.eps
        Hmat[1,0] = self.Delta
        Hmat[0,1] = self.Delta
        return Hmat

    def pos_matrix(self): #perturbation matrix
        s_mat = np.zeros((self.ns,self.ns),dtype=complex)
        s_mat[0,0] = 1
        s_mat[1,1] = 1
        return s_mat

    # Generates the initial system density matrix and the partition function
    def initcond(self,beta): #(here we just put the system in the excited state)
        rho_s0 = np.zeros((self.ns,self.ns),dtype=complex)
        rho_s0[1,1] = 1
        Zs = 1
        return rho_s0,Zs

    # Generate the correlation function (this is hardcoded)
    def corr(self,rho_s,t): # (return population of the excited state)
        return rho_s[1,1],t

    def analytic_uncoupled(self,beta,t_arr=np.arange(0,10,0.1)):# calculatees the analytic solution for the uncoupled system
        H= self.H_matrix()
        vals,vecs = np.linalg.eigh(H)
        phi_0=np.zeros((self.ns),dtype=complex)
        phi_0[1]=1
        phi_0 = vecs.T@phi_0 #transform into eigenbasis
        pop = np.zeros_like(t_arr)
        for i,t in enumerate(t_arr):
            K = np.diag(np.exp(-1j*vals*t))
            phi_t = K@phi_0
            pop[i] = np.abs((vecs@phi_t)[1])**2
        return t_arr,pop