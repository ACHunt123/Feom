import numpy as np
import sys  
from scipy import special
import matplotlib.pyplot as plt


### Generates the Hamiltonian matrix spin boson model in diabatic representation
class Spin_boson:
    def __init__(self,Delta,eps):
        self.m = 1
        self.ns = 2 # Hardcoded as spinboson is twolevel
        self.Delta = Delta
        self.eps = eps



    # Generates the Hamiltonian matrix in its eigenbasis
    def H_matrix(self):
        Hmat = np.zeros((self.ns,self.ns),dtype=complex)
        Hmat[0,0] = -self.eps
        Hmat[1,1] = self.eps
        Hmat[1,0] = self.Delta
        Hmat[0,1] = self.Delta
        return Hmat

    def pos_matrix(self):
        psi_ns, E_ns, x_arr = self.eigenstates()
        # Generate the position matrix in energy basis
        pos_matrix = np.zeros((self.ns,self.ns),dtype=complex)                   
        for m in range(0,self.ns):
            for n in range(0,self.ns):
                pos_matrix[m,n] = np.sum(x_arr[:]*np.conj(psi_ns[:,m])*psi_ns[:,n]*self.dx)
                # pos_matrix[m,n] = np.vdot(np.conj(psi_ns[:,m]),x_arr[:]*psi_ns[:,n])*self.dx
        return pos_matrix#, x_arr, dx

    # Generates the initial system density matrix and the partition function
    def initcond(self,beta): #(here we just put the system in the excited state)
        rho_s0 = np.zeros((self.ns,self.ns),dtype=complex)
        rho_s0[1,1] = 1
        Zs = 1
        return rho_s0,Zs

    # Generate the correlation function (this is hardcoded)
    def corr(self,rho_s,t): # (return population of the excited state)
        return rho_s[1,1],t