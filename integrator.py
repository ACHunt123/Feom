import numpy as np
import scipy, sys
import matplotlib.pyplot as plt
from hashmap import generateHashmap
#
#   This is the integrator for the HEOM
#
#   it is done in the energy eigenbasis
class Integrator:
    def __init__(self,gam_ks,C_ks,ns,Imax,H_mat,s_mat,K,hbar,L,lowTcoef):
        self.ns = ns
        self.Imax = Imax
        self.L = L
        self.K = K
        self.hbar = hbar
        self.lowTcoef=lowTcoef
        # The formatting of the Hashmaps are as follows:
        # I2ind[I] : int I -> list of ints corresponding to the BCF indecies
        # ind2I[ind] :tuple of ints -> int I
        # This is done because lists are not hashable, and tuples are   
        self.I2ind , self.ind2I = generateHashmap(K,L)
        self.H_mat = H_mat
        self.s_mat = s_mat # perturbation operator (could be q for example) 
        # Bath coefficients
        self.gam_ks = gam_ks
        self.C_ks = C_ks
        # objects in RK4 calculation
        self.k1 = np.zeros((ns,ns,Imax),dtype=complex)
        self.k2 = np.zeros((ns,ns,Imax),dtype=complex)
        self.k3 = np.zeros((ns,ns,Imax),dtype=complex)
        self.k4 = np.zeros((ns,ns,Imax),dtype=complex)

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
        elif self.tier(ind) > self.L:
            return -1
        # Otherwise, the function find the index of the ADO with the same elements as ind, but with the kth element changed to nkpm and return it
        else:
            ind = str(ind).replace('[','').replace(']','').replace(' ','')
            return self.ind2I[ind]  

    def commutator(self,A,B):
        return (A@B - B@A)

    # gives the unperturbed Liouvillian acting on the ADO = i/hbar [H,rho]
    def L0(self,rho):
        return 1.j/self.hbar * self.commutator(self.H_mat,rho)

    # gives the gradient of the density matrix with respect to time
    def drhodt(self,rho, gradient):
        gradient[:,:,:] = 0

        for I in range(0,self.Imax):
            n_ks = np.array(self.I2ind[I]) # the indexes of the ADO - will be used as coefficients 

            gradient[:,:,I] = -self.L0(rho[:,:,I]) - np.sum(n_ks*self.gam_ks)*rho[:,:,I]

            if self.lowTcoef!=0 : # if the bath is in the matsubara mode employ the Ishizaki-Tanimura method
                gradient[:,:,I] -= self.commutator(self.s_mat,self.commutator(self.s_mat,rho[:,:,I])) * self.lowTcoef

            for k in range(self.K+1): # as there are K+1 terms in the BCF

                # the ADOs that are one element different from the current ADO
                I_nkp1 = self.I_nk_plusminus(I,k,+1)
                I_nkm1 = self.I_nk_plusminus(I,k,-1)

                if(0):
                    print('I:',I,'===',self.I2ind[I], '-- k =',k)
                    print('I_nkp1:',self.I2ind[I_nkp1]) if I_nkp1 != -1 else print('I_nkp1:',I_nkp1,'===','not here')
                    print('I_nkm1:',self.I2ind[I_nkm1]) if I_nkm1 != -1 else print('I_nkm1:',I_nkm1,'===','not here')
                    print('---')

                Ck = self.C_ks[k]
                nk = n_ks[k]
                absCk = np.abs(Ck)

                # The gradient for the +1 terms [may not exist]
                if I_nkp1 != -1:
                    gradient[:,:,I] -= 1.j/self.hbar *np.sqrt((nk+1)*absCk) * self.commutator(self.s_mat,rho[:,:,I_nkp1])

                # The gradient for the -1 terms [may not exist]
                if I_nkm1 != -1:
                    gradient[:,:,I] -= 1.j/self.hbar * np.sqrt(nk/absCk) * (Ck*self.s_mat@rho[:,:,I_nkm1] - np.conj(Ck)*rho[:,:,I_nkm1]@self.s_mat)

    # RK4 step
    def rk4_step(self,x0,dt):
        self.drhodt(x0,gradient=self.k1)
        self.k1 = self.k1 * (dt/2)

        self.drhodt(x0+self.k1,gradient=self.k2)
        self.k2 = self.k2 * (dt/2)

        self.drhodt(x0+self.k2,gradient=self.k3)
        self.k3 = self.k3 * (dt)
        
        self.drhodt(x0+self.k3,gradient=self.k4)
        return x0 + dt/6.*self.k4 + (2./3.)*self.k2 + (self.k3 + self.k1)/3.

