import numpy as np
import sys, os
import matplotlib.pyplot as plt
from hashmap import generateHashmap,Convert_to_list
import Heom.fort.executables.propagator as prop
npF = np.asfortranarray # Aliasing to make the code more legible

#
#   This is the integrator for the HEOM - done with scaled matrices st dx = 1 (or discrete basis)
#
class Integrator:
    def __init__(self,gam_ks,C_ks,ns,Imax,H_mat,s_mat,K,hbar,L,lowTcoef,N_nonmats):
        self.ns = ns
        self.Imax = Imax
        self.L = L
        self.K = K
        self.hbar = hbar
        self.lowTcoef=lowTcoef
        self.N_nonmats = N_nonmats
        # The formatting of the Hashmaps are as follows:
        # I2ind[I] : int I -> list of ints corresponding to the BCF indecies
        # ind2I[ind] :tuple of ints -> int I
        # This is done because lists are not hashable, and tuples are   
        self.I2ind , self.ind2I = generateHashmap(K,L,N_nonmats) #hash map from the index of the ADO to the index of the BCF
        self.ADO_index, self.I0s = Convert_to_list(self.I2ind) #new indexing 
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

    ### Slower but more verbose version of function that can be put into FORTRAN
    # def I_nk_plusminus(self,I,k,pm):
    #     index = self.ADO_index[I,:].copy() 
    #     tier0 = self.tier(index) # Initial tier of the ADO
    #     # Calculate the new index
    #     index[k] += pm
    #     # Check if the new index is valid
    #     if  self.tier(index)> self.L :
    #         return -1
    #     if index[k]<0 :
    #         return -1
    #     # Find the place to search for the ADO
    #     # I0/I1 = minimum/minimum index of where to look for the next ADO
    #     if pm == +1:
    #         I0 = self.I0s[tier0+1]
    #         I1 = self.I0s[tier0+2]
    #     else:
    #         I0 = self.I0s[tier0-1]
    #         I1 = self.I0s[tier0]

    #     def find_index(indices, index):
    #         for i in range(len(indices)):
    #             if np.all(indices[i] == index):
    #                 return i
    #         raise ValueError('Index not found')

    #     return find_index(self.ADO_index[I0:I1,:],index)+I0

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
    def rk4_step(self,x0,dt,FORTRAN=0):
        if not FORTRAN: # RK$ using python
            self.drhodt(x0,gradient=self.k1)
            self.k1 = self.k1 * (dt/2)

            self.drhodt(x0+self.k1,gradient=self.k2)
            self.k2 = self.k2 * (dt/2)

            self.drhodt(x0+self.k2,gradient=self.k3)
            self.k3 = self.k3 * (dt)
            
            self.drhodt(x0+self.k3,gradient=self.k4)
            return x0 + dt/6.*self.k4 + (2./3.)*self.k2 + (self.k3 + self.k1)/3.
        else: # Propagation using FORTRAN
            # Format all of the data
            x0fort =np.zeros((self.Imax,self.ns,self.ns),dtype=complex,order='F')
            for I in range(self.Imax):
                x0fort[I,:,:] = npF(x0[:,:,I])
            # Propagate the density matrix
            prop.prop_subroutines.vvstep(x0fort,npF(self.ADO_index),npF(self.I0s),npF(self.gam_ks),
                npF(self.C_ks),npF(self.H_mat),npF(self.s_mat),self.K,self.hbar,dt,
                self.lowTcoef,self.N_nonmats,imax=self.Imax,l=self.L,ns=self.ns)
            # Reformat the density matrix for python
            for I in range(self.Imax):
                x0[:,:,I] = x0fort[I,:,:]
            return x0

    def generate_input_files():
        os.system('mkdir tmp') #make a temporary directory to store the input files
        # Format all of the data
        x0fort =np.zeros((self.Imax,self.ns,self.ns),dtype=complex,order='F')
        for I in range(self.Imax):
            x0fort[I,:,:] = npF(x0[:,:,I])
        ADO_index = npF(self.ADO_index)
        I0s = npF(self.I0s)
        gam_ks = npF(self.gam_ks)
        C_ks = npF(self.C_ks)
        H_mat = npF(self.H_mat)
        s_mat = npF(self.s_mat)
        K = self.K
        L = self.L
        hbar = self.hbar
        lowTcoef = self.lowTcoef
        N_nonmats = self.N_nonmats
        Imax = self.Imax
        ns = self.ns
        return
