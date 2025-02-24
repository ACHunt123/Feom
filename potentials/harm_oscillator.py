import numpy as np
import sys  
from scipy import special
import matplotlib.pyplot as plt


### Generates the Hamiltonian matrix for a harmonic oscillator in the position basis
class Harmonic_oscillator:
    def __init__(self,ns=10):
        ### Hamiltonian and system setup - harmonic oscillator - same as Adam's code
        m =1741.1
        d0 = 0.18748
        alpha = 1.1605
        diff = 2 * d0 *  alpha**2
        const = d0 *  alpha**2
        omega = (diff/m)**(0.5)

        ### Position basis parameters for calculation of matrices
        xmin = -5
        xmax = 5
        dx = 0.01
        ns = 10         # number of eigenstates to be propagated


        self.m = m
        self.omega = omega
        self.hbar = 1
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
            E_ns[i] = self.hbar*self.omega*(i+0.5)    # No lamb shift here

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
        return pos_matrix

    # Generates the initial system density matrix and the partition function
    def initcond(self,beta):
        # Get eigenstates (alternatively we could use DVRs)
        psi_ns, E_ns, x_arr = self.eigenstates()
        delEs = E_ns - E_ns[0] # the ground state is the zero of energy
        rho_s = np.diag(np.exp(-beta*delEs))
        Zs = np.trace(rho_s)
        rho_s0 = rho_s@self.pos_matrix()/Zs
        return rho_s0,Zs

    # Generate the correlation function (this is hardcoded)
    def corr(self,rho_s,t):
        return np.trace(rho_s[:,:]@self.pos_matrix()),t

    
    def analytic_uncoupled(self,beta,t_arr=np.arange(0,10,0.1)):# calculatees the analytic solution for the uncoupled system
        x = np.exp(beta*self.hbar*self.omega/2)
        xm1 = x**-1
        Css_analyt_re = ((self.hbar/(2*self.m*self.omega)) * (x+xm1)/(x-xm1))*np.cos(self.omega*t_arr)

        return t_arr,Css_analyt_re