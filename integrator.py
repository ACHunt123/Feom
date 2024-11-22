import numpy as np
import scipy
import matplotlib.pyplot as plt
from hashmap import generateHashmap
#
#   This is the integrator for the HEOM
#
class Integrator:
    def __init__(self,ns,ds,Imax,H_mat,s_mat,K,hbar,L):
        self.ns = ns
        self.Imax = Imax
        self.L = L
        self.K = K
        self.hbar = hbar
        # The formatting of the Hashmaps are as follows:
        # I2ind[I] : int I -> list of ints corresponding to the BCF indecies
        # ind2I[ind] :tuple of ints -> int I
        # This is done because lists are not hashable, and tuples are   
        self.I2ind , self.ind2I = generateHashmap(K,max_N)
        self.H_mat = H_mat
        self.s_mat = s_mat
        self.ds = ds

    # returns the tier of the ADO from its list of indices
    def tier(self,ind):
        return np.sum(ind)

    # returns the index of the ADO that is one element different from the input ADO
    def I_nk_plusminus(self,I,k,pm):
        index = self.I2ind[I] #copy the list of indexes of the ADO
        ind = index.copy() #copy the list of indexes of the ADO

        nkpm = ind[k]+pm
        ind[k] = nkpm   #update the index of the ADO with the kth element changed to nkpm

        # if nkpm is negative, the ADO does not exist the function returns -1
        if nkpm < 0:
            return -1
        # if the tier of the ADO index does not exist the function returns -1
        elif self.tier(ind) >= self.L:
            return -1
        # Otherwise, the function find the index of the ADO with the same elements as ind, but with the kth element changed to nkpm and return it
        else:
            ind = str(ind).replace('[','').replace(']','').replace(' ','')
            return self.ind2I[ind]  

    def commutator(self,A,B):
        return (A@B - B@A)*self.ds

    # gives the unperturbed Liouvillian acting on the ADO = i/hbar [H,rho]
    def L0(self,rho):
        return 1.j/self.hbar * self.commutator(self.H_mat,rho)

    # gives the gradient of the density matrix with respect to time
    def drhodt(self,rho):
        gradient = np.zeros((self.ns,self.ns,self.Imax),dtype=complex)

        for I in range(0,self.Imax):
            # n_ks = np.array(self.I2ind[I]) # the indexes of the ADO - will be used as coefficients 

            gradient[:,:,I] = -self.L0(rho[:,:,I]) #+ np.sum(n_ks*gam_ks)*rho[:,:,I]


            # for k in range(K):
            #     # the ADOs that are one element different from the current ADO
            #     I_nkp1 = self.I_nk_plusminus(I,k,+1)
            #     I_nkm1 = self.I_nk_plusminus(I,k,-1)

            #     # The gradient for the +1 terms [may not exist]
            #     if I_nkp1 != -1:
            #         gradient[:,:,I] += 1.j/self.hbar * self.commutator(self.s_mat,rho[:,:,I_nkp1]) * np.sqrt(self.C_ks[k]*(self.n_ks[k]+1)) ###redo#

            #     # The gradient for the -1 terms [may not exist]
            #     if I_nkm1 != -1:
            #         gradient[:,:,I] += 1.j/hbar * self.commutator(self.s_mat,rho[:,:,I_nkp1]) * np.sqrt(C_ks[k]*(n_ks[k]+1)) ###redo#

        return gradient

    # RK4 step
    def rk4_step(self,x0,dt):
        k1 = self.drhodt(x0)
        k2 = self.drhodt(x0 + 0.5*dt*k1)
        k3 = self.drhodt(x0 + 0.5*dt*k2)
        k4 = self.drhodt(x0 + dt*k3)
        return x0 + dt/6*(k1 + 2*k2 + 2*k3 + k4)

if(0): #old rk4 testing code
    p0 = np.ones(2)
    A = np.array([[0,1],[-1,2]])
    def grad(p):
        return A.T@p


    dt = 0.1
    p = p0
    ps = np.zeros((100,2))
    analyt_ps = np.zeros((100,2))
    t = []
    for i in range(100):
        ti=i*dt
        p_dp = rk4_step(p,dt)
        ps[i,:] = p
        t.append(ti)
        analyt_ps[i,:] = p0@scipy.linalg.expm(A*ti)
        p = p_dp


    plt.plot(t,ps[:,0])
    plt.plot(t,ps[:,1])
    plt.plot(t,analyt_ps[:,0],'--')
    plt.plot(t,analyt_ps[:,1],'--')
    plt.show()
