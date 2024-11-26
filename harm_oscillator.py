import numpy as np
import sys  
from scipy import special
import matplotlib.pyplot as plt


### Generates the Hamiltonian matrix for a harmonic oscillator in the position basis
class Harmonic_oscillator:
    def __init__(self,m,omega,hbar,xmin,xmax,dx,ns):
        self.m = m
        self.omega = omega
        self.hbar = hbar
        self.xmin = xmin
        self.xmax = xmax
        self.dx = dx
        self.ns = ns
        self.nx = int((self.xmax-self.xmin)/self.dx) # number of x points
        

        
    # Generates the eigenstates of the Hamiltonian in the posision basis - to be replaced with DVRs
    def eigenstates(self):
        x_arr= np.zeros(self.nx)      # x points
        for i in range(self.nx):
        # + 1 is there for consitency with adams code
            x_arr[i] = self.dx*(i + 1 + int(self.xmin/self.dx))

        psi_ns = np.zeros((self.nx,self.ns),dtype=complex)  # the eigenvectors of the Hamiltonian
        E_ns = np.zeros(self.ns)                       # the eigenvalues of the Hamiltonian
        for i in range(0,self.ns):
            alph = self.m*self.omega/self.hbar
            prefac = 1/np.sqrt(2.**i * np.math.factorial(i))   * (alph/np.pi)**0.25
            Hi = special.hermite(i)
            psi_ns[:,i] = prefac * Hi(np.sqrt(alph)*x_arr) * np.exp(-alph*x_arr**2/2)
            E_ns[i] = self.hbar*self.omega*(i)    # No lamb shift here

        if(0):# plot the eigenstates
            fig = plt.figure()
            ax = fig.add_subplot(111)
            for i in range(0,self.ns):
                ax.plot(x_arr,psi_ns[:,i],label='En='+str(E_ns[i]))
            plt.legend()
            plt.show()

        if(0):# Test that the eigenstates are orthonormal
            for i in range(0,self.ns):
                for j in range(0,self.ns):
                    print(np.round(np.sum(psi_ns[:,i]*np.conj(psi_ns[:,j])*dx),4),i,j)
        
        return psi_ns, E_ns, x_arr

    # Generates the Hamiltonian matrix in its eigenbasis
    def H_matrix(self):
        psi_ns, E_ns, x_arr = self.eigenstates()
        return np.diag(E_ns)

    def pos_matrix(self):
        psi_ns, E_ns, x_arr = self.eigenstates()
        # Generate the position matrix in energy basis
        pos_matrix = np.zeros((self.ns,self.ns),dtype=complex)                   
        for m in range(0,self.ns):
            for n in range(0,self.ns):
                pos_matrix[m,n] = np.sum(x_arr[:]*np.conj(psi_ns[:,m])*psi_ns[:,n]*self.dx)
        return pos_matrix#, x_arr, dx

    # Generates the initial system density matrix in the position basis
    def rho0(self,beta):
        # Get eigenstates (alternatively we could use DVRs)
        psi_ns, E_ns, x_arr = self.eigenstates()
        return np.diag(np.exp(-beta*E_ns))