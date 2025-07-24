import numpy as np
import sys  
from scipy import special
import matplotlib.pyplot as plt


### Generates the Hamiltonian matrix for a harmonic oscillator in the position basis
class Harmonic_oscillator:
    def __init__(self,params):
        self.m=params.m
        self.omega=params.omega
        self.beta=params.beta
        self.hbar=params.hbar
        self.xmin=params.xmin
        self.xmax=params.xmax
        self.dx=params.dx
        self.ns=params.ns
        self.nx = int((self.xmax-self.xmin)/self.dx) # number of x points
        ### Generate the Hamiltonian matrix in the eigenbasis
        self.H_matrix()
        ### Generate the position matrix in the eigenbasis
        self.s_mat = self.pos_matrix()
        # plot the initial conditions matrix
        if(0):
            fig = plt.figure()
            ax = fig.add_subplot(111)
            rho,z= self.initcond()
            print(z)
            print(rho)
            ax.imshow(np.abs(rho), interpolation='nearest')
            plt.title('Initial Conditions Matrix')
            plt.xlabel('Energy Basis')
            plt.ylabel('Energy Basis')
            plt.show()
            sys.exit()
        return

        
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
                    print(np.round(np.sum(psi_ns[:,i]*np.conj(psi_ns[:,j])*self.dx),4),i,j)
        
        return psi_ns, E_ns, x_arr

    # Generates the Hamiltonian matrix in its eigenbasis
    def H_matrix(self):
        psi_ns, E_ns, x_arr = self.eigenstates()
        self.H_mat = np.diag(E_ns)
        return 

    def pos_matrix(self):
        psi_ns, E_ns, x_arr = self.eigenstates()
        # Generate the position matrix in energy basis
        pos_matrix = np.zeros((self.ns,self.ns),dtype=complex)                   
        for m in range(0,self.ns):
            for n in range(0,self.ns):
                pos_matrix[m,n] = np.sum(x_arr[:]*np.conj(psi_ns[:,m])*psi_ns[:,n]*self.dx)
        return pos_matrix

    # Generates the initial system density matrix and the partition function
    def initcond(self):
        # Get eigenstates (alternatively we could use DVRs)
        psi_ns, E_ns, x_arr = self.eigenstates()
        delEs = E_ns - E_ns[0] # the ground state is the zero of energy
        rho_s = np.diag(np.exp(-self.beta*delEs))
        Zs = np.trace(rho_s)
        rho_s0 = rho_s@self.pos_matrix()/Zs # the intial condition is rho q 
        # rho_s0 = rho_s/Zs # the intial condition is rho  
        return rho_s0,Zs

    # Generate the correlation function (this is hardcoded)
    def corr(self,rho_s,t):
        return np.trace(rho_s[:,:]@self.pos_matrix()),t

    # Calculate the analytic solution for the uncoupled system
    def analytic_uncoupled(self,t_arr=np.arange(0,10,0.1)):
        x = np.exp(self.beta*self.hbar*self.omega/2)
        xm1 = x**-1
        Css_analyt_re = ((self.hbar/(2*self.m*self.omega)) * (x+xm1)/(x-xm1))*np.cos(self.omega*t_arr)
        return t_arr,Css_analyt_re
    
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
        # calculate <q> 
        processed_data = np.zeros((len(t),5),dtype=complex)

        processed_data[:,0] = t
        processed_data[:,1] = np.einsum('inm,mn->i', rho,self.s_mat)
        for it in range(len(t)):
            processed_data[it,2] = np.min(np.diagonal(rho[it,:,:]))  # min diagonal element
            processed_data[it,3] = np.sum(np.diagonal(rho[it,:,:]))  # trace of the density matrix

        data_labels = '\n,Time /a.u. <q> min.diag. trace\n'
        return processed_data, initial_header+data_labels