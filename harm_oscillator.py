import numpy as np
from scipy import special
import matplotlib.pyplot as plt


### Generates the Hamiltonian matrix for a harmonic oscillator in the position basis

# Generates the eigenstates of the Hamiltonian
def eigenstates(ns, m, omega, hbar, nx, xmin, xmax):
    nvecs = 25
    x_arr= np.linspace(xmin,xmax,nx)
    dx = x_arr[1]-x_arr[0]

    psi_ns = np.zeros((ns,nvecs),dtype=complex)  # the eigenvectors of the Hamiltonian
    E_ns = np.zeros(nvecs)                       # the eigenvalues of the Hamiltonian
    for i in range(0,nvecs):
        alph = m*omega/hbar
        prefac = 1/np.sqrt(2.**i * np.math.factorial(i))   * (alph/np.pi)**0.25
        Hi = special.hermite(i)
        psi_ns[:,i] = prefac * Hi(np.sqrt(alph)*x_arr) * np.exp(-alph*x_arr**2/2)
        E_ns[i] = hbar*omega*(i+0.5)

    if(0):# plot the eigenstates
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for i in range(0,25):
            ax.plot(x_arr,psi_ns[:,i],label='En='+str(E_ns[i]))
        plt.legend()
        plt.show()
    
    return psi_ns, E_ns, x_arr, dx

# Generates the Hamiltonian matrix in the position basis
def H_matrix(ns=150, m=1, omega=1, hbar=1, nx=150, xmin=-14, xmax=14):
    # Get eigenstates (alternatively we could use DVRs)
    psi_ns, E_ns, x_arr, dx = eigenstates(ns, m, omega, hbar, nx, xmin, xmax)
    # Generate the Hamiltonian matrix in the position basis
    H = np.zeros((ns,ns),dtype=complex)                   
    for s in range(0,ns):
        for sp in range(0,ns):
            H[s,sp] = np.sum(E_ns*psi_ns[s,:]*np.conj(psi_ns[sp,:]))
    
    if(0):# plot a surface plot of the Hamiltonian
        fig = plt.figure()
        ax = fig.add_subplot(111)
        X, Y = np.meshgrid(x_arr, x_arr)
        ax.contourf(X, Y, H, cmap='viridis')
        plt.show()
    
    return H, x_arr, dx

# Generates the initial system density matrix in the position basis
def rho0(beta=6, ns=150, m=1, omega=1, hbar=1, nx=150, xmin=-14, xmax=14):
    # Get eigenstates (alternatively we could use DVRs)
    psi_ns, E_ns, x_arr, dx = eigenstates(ns, m, omega, hbar, nx, xmin, xmax)
    # Generate the Hamiltonian matrix in the position basis
    rho = np.zeros((ns,ns),dtype=complex)                   
    for s in range(0,ns):
        for sp in range(0,ns):
            rho[s,sp] = np.sum(np.exp(-beta*E_ns)*psi_ns[s,:]*np.conj(psi_ns[sp,:]))

    return rho