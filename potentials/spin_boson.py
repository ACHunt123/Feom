import numpy as np
import sys  
from scipy import special
import matplotlib.pyplot as plt


### Generates the Hamiltonian matrix spin boson model in diabatic representation
# Note adiabatic repn is without _'s and diabatic repn is with _'s on matrices
class Spin_boson:
    def __init__(self,params):
        assert params.ns == 2, "Spin boson model is hardcoded for two levels"
        self.ns = params.ns
        self.Delta = params.Delta
        self.eps = params.eps
        ### Generate the Hamiltonian matrix in the diabatic basis
        self.H_matrix()
        ### compute the diabatic<->adiabatic transformation matrix ONCE
        self.eigs,self.Uda = np.linalg.eigh(self.H_mat) # diagonalize the Hamiltonian
        ### generate the perturbation matrix in the diabatic basis
        self.s_mat = self.pos_matrix()
        return

    # Generates the Hamiltonian matrix in the diabatic basis
    def H_matrix(self):
        self.H_mat = np.zeros((self.ns,self.ns),dtype=complex)
        self.H_mat[0,0] = -self.eps
        self.H_mat[1,1] = self.eps
        self.H_mat[1,0] = self.Delta
        self.H_mat[0,1] = self.Delta
        return 

    def pos_matrix(self): # perturbation matrix in diabatic basis
        q_mat = np.zeros((self.ns,self.ns),dtype=complex)
        q_mat[0,0] = -1
        q_mat[1,1] = 1
        return q_mat

    # Generates the initial system density matrix and the partition function
    def initcond(self): #(here we just put the system in the excited state)
        rho_s0 = np.zeros((self.ns,self.ns),dtype=complex)
        rho_s0[1,1] = 1
        Zs = 1
        return rho_s0,Zs

    # Generate the correlation function (this is hardcoded)
    def corr(self,rho_s,t): # (return population of the excited diabatic state)
        return rho_s[1,1],t

    # Calculatees the analytic solution for the uncoupled system
    def analytic_uncoupled(self,t_arr=np.arange(0,10,0.1)):
        phi_0=np.zeros((self.ns),dtype=complex)
        phi_0[1]=1
        phi_0 = self.Uda.T@phi_0 #transform into eigenbasis
        pop = np.zeros_like(t_arr)
        for i,t in enumerate(t_arr):
            K = np.diag(np.exp(-1j*self.eigs*t))
            phi_t = K@phi_0
            pop[i] = np.abs((self.Uda@phi_t)[1])**2
        return t_arr,pop
       
    # Format FORTRAN output and calculate common observables
    def format_output(self,data,initial_header):
        # Collect the data
        t= data[:,0]
        re_rho = data[:,1:1+self.ns**2]
        im_rho = data[:,1+self.ns**2:]
        rho = re_rho + 1.j*im_rho
        rho = rho.reshape((len(t),self.ns,self.ns), order='F')  # reshape the data to be a 3D array
        # format the data to calculate <s_z>, <s_y> and <s_x> and others
        t = data[:,0]
        rho11 = rho[:,1,1]  
        rho00 = rho[:,0,0]
        rho10 = rho[:,1,0]
        rho01 = rho[:,0,1]

        processed_data = np.zeros((len(t),5),dtype=complex)
        processed_data[:,0] = t
        processed_data[:,1] = rho11 - rho00  # <s_z>
        processed_data[:,2] = 1.j*(rho10 - rho01)  # <s_y>
        processed_data[:,3] = (rho10 + rho01)  # <s_x>   
        processed_data[:,4] = (1+(rho11 - rho00))/2  # Site 1 population
        data_labels = '\n,Time /a.u. <s_z> <s_y> <s_x> (1+<s_z>)/2'
        return processed_data, initial_header+data_labels
