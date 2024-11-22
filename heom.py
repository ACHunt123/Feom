from hashmap import  total_length
from harm_oscillator import H_matrix, rho0
from debye_bath import Debye_bath
from integrator import Integrator
import numpy as np
import scipy
import sys
import matplotlib.pyplot as plt
#
# Basic HEOM code
#


### Parameters
beta = 10
hbar=1
ns = 33         # number of states/xtics (same thing really as we are in the position basis)
K = 1           #  the number of elements in the BCFs
L = 1           # the depth of the ADO expansion

### Bath parameters - Debye bath
bathmode = 'matsubara'
eta=0.1
gam=1


### Get bath coefficients
bath = Debye_bath(eta,gam,beta,hbar,K,bathmode)
C_ks,gam_ks = bath.get_coeffs()


### Hamiltonian and system setup - harmonic oscillator
m =1
omega = 1

nx = ns          #number of xtics = number of states
xmin = -10
xmax = 10


H_mat, x_arr, ds = H_matrix(ns,m,omega,hbar,nx,xmin,xmax)  # Hamiltonian matrix

s_mat = np.diag(x_arr) 
rho_s0 = rho0(beta,ns,m,omega,hbar,nx,xmin,xmax)           # initial system density matrix = e^(-beta*H_s) 

#Trace operator
def trace(matrix): return np.sum(np.diag(matrix))*ds    # the extra factor of ds is because we are in the position basis and it is continuous

#partition function
Zs = trace(rho_s0)                                    # partition function - calculated without lamb shift



if(1):  # tests
    # Test that the partition function is correct
    theta=hbar*omega*beta
    print(1/(1-np.exp(-theta)),'analytic partition function')
    print(Zs,'partition function calculated','\n')

    # Test that the thermal expectation value of q^2 is correct
    x = np.exp(-beta*hbar*omega)
    print((hbar/(2*m*omega)) * (1+x)/(1-x), 'analytic <q^2>')

    # Calculate the thermal expectation value of q^2 for the system
    print(trace(rho_s0@s_mat**2*ds)/Zs,'<q^2> calculated')

    print(trace(rho_s0@s_mat**2)/Zs,'<q^2> calculated WITHOUT EXTRA FACTOR OF ds') #this is incorrect but seems to work


    # Test that the thermal density matrix is stationary 
    # print(((rho_s0@H_mat - H_mat@rho_s0)*ds).max(),'maximal deviation from stationarity')

sys.exit()

### Variables and arrays
Imax = total_length(K,L)                # the total number of ADOs
# Iprop = #numbers of ADOs to propagete on the current time step

rho = np.zeros((ns,ns,Imax),dtype=complex)  # holds all of the ADOs. rho[s,s',0] is the system density matrix
rho[:,:,0] = rho_s0@s_mat*ds                # initial condition = e^-beta*H_s * s

#setup the integrator
integrator = Integrator(ns,ds,Imax,H_mat,s_mat,K,hbar,L)

### Propagation
tmax = 10
dt=0.01
nt =int(tmax/dt)
t_arr = np.linspace(0,tmax,nt)
dt = t_arr[1]-t_arr[0]

### Plotting
fig, (ax1, ax2) = plt.subplots(1,2)
plt.ion()

# plot the analytical solution for uncoupled harmonic oscillator
x = np.exp(beta*hbar*omega/2)
xm1 = x**-1
Css_analyt_re = ((hbar/(2*m*omega)) * (x+xm1)/(x-xm1))*np.cos(omega*t_arr)
ax2.plot(t_arr,Css_analyt_re,'--',color='red',label='analytical')
ax2.set_ylim(-1,1)
Css = np.zeros_like(t_arr)
line, = ax2.plot(t_arr,np.zeros_like(t_arr),'b',label='numerical')

for i in range (0,1000):

    rho = integrator.rk4_step(rho,dt)

    # print(rho[:,:,0]@x_*ds)
    Css[i] = trace(rho[:,:,0]@s_mat*ds)/Zs #there is a factor of ds on the bottom missing
    # plot the density matrix
    if(i%10==0):
        line.remove()
        ax1.clear()
        line, = ax2.plot(t_arr[:i],Css[:i],'b',label='numerical')
        ax1.contourf(x_arr,x_arr,np.real(rho[:,:,0]),cmap='viridis')
        # print(trace(rho[:,:,0].T@np.diag(x_arr**2)*ds)/np.trace(rho[:,:,0]*ds))
        plt.pause(0.01)
plt.show()

### Functions


