import numpy as np
import matplotlib.pyplot as plt
import sys

# A class for the Debye bath
''' A class to represent the Debye bath for the FEOM code.

eta : Coupling strength
gam : Cuttoff frequency

The spectral density is given by:
J(w) = \frac{\eta\gamma\omega}{\omega^2 + \gamma^2}
using the discretized definition, giving coefficients:
J(w) = (\pi/2) * \sum_{\alpha} \frac{c_\alpha^2}{m_\alpha \omega_\alpha} \delta(w - w_\alpha)
'''


class Debye_bath():
    def __init__(self,params):
        # Bathmode and settings
        self.bathmode = params.bathmode
        self.L = params.L                      # max tier of the ADOs
        # General parameters
        self.eta = params.eta
        self.gam = params.gam
        self.beta = params.beta
        self.hbar = params.hbar
        # Paramaters for bath indexing
        self.N_nonmats = 1                          # number of exponential modes in BCF that are NOT matsubara terms (the Temp. ind. exp.)
        self.mu = params.K                          # number of pairs of matsubara modes/ r.p. modes, each pair gives a single exponential term
        self.N_exp = self.N_nonmats + self.mu       # number of exponential terms in the BCF [Temp ind. Exponential, <--- Matsubara Exponentials --->]
        self.N_exp_prop = self.N_nonmats + self.mu  # number of exponentials that we are propogating (may be different to the number used in the coth decomposition)
        self.N_mds = 2*self.mu+1                    # number of individual beads or matsubara modes (ODD)
        #
        self.mode= params.bathmode
        self.C0hot = self.eta/self.beta -1.j*self.hbar*self.eta*self.gam/2 # C_0 with no matsubara terms
        ### Calculate the C_ks and gam_ks for the bath and add to the class
        self.get_coeffs()
        ### Calculate the coefficients C_U, c_D_LEFT, c_D_RIGHT for the bath (that are used in the FEOM code)
        self.get_C_UDs()

    def J(self,w,plotme=False,ax=plt):
        w = np.linspace(0,2,1000)
        Jw = self.eta*self.gam*w/(w**2+self.gam**2)
        if plotme : 
            ax.plot(w,Jw)
        return Jw,w

    # Calculate the C_ks and gam_ks for a given set of ws
    def calc_coefs(self,ws):
        d = np.zeros(self.N_exp)
        ### Calculate the 0th (non-matsubara) term
        if self.mode in ['nbead','nmats']:              # give the finite mode prefactor
            d0sum = np.sum(1/(self.gam**2*np.ones_like(ws[1:]) - ws[1:]**2))
            d[0] = (self.hbar*self.eta*self.gam/2) * (2/(self.beta*self.hbar*self.gam) + (4*self.gam/(self.beta*self.hbar))*d0sum)
        elif self.mode == 'matsubara': # give the infinite mode prefactor
            d[0] = (self.hbar*self.eta*self.gam/2) /np.tan(self.beta*self.hbar*self.gam/2)
        else:
            raise ValueError('Invalid mode')
        ### Calculate the rest of the terms (matsubara terms)
        d[1:] = -(2*self.eta*self.gam/self.beta) * ws[1:]/(self.gam**2*np.ones_like(ws[1:]) - ws[1:]**2)
        dI = -self.hbar*self.eta*self.gam/2

        # Calculate the C_ks
        C_ks = np.zeros(self.N_exp,dtype=complex)
        C_ks[0] = (d[0] + dI*1.j)
        C_ks[1:] = d[1:]

        #Calculate the gam_ks
        gam_ks = np.zeros(self.N_exp,dtype=complex)
        gam_ks[0] = self.gam
        gam_ks[1:] = ws[1:]
        
        # Calculate the low temperature coefficient
        if self.mode in ['nmats','nbead']: 
            self.lowTcoef = 0
        elif self.mode == 'matsubara':  # The Ishizaki-Tanimura terminator coefficient
            self.lowTcoef = self.eta* ((1/(2*self.hbar))* ((1/(np.tan(self.beta*self.hbar*self.gam/2))) - (2/(self.beta*self.hbar*self.gam))))  ### Terms without removing of the matsubara terms that have been included
            self.lowTcoef = self.lowTcoef -  self.eta*(2*self.gam/(self.beta*self.hbar**2))*np.sum(1/(self.gam**2*np.ones_like(ws[1:]) - ws[1:]**2))  ### remove the Matsubara terms that have been explicitly included

        self.C_ks = C_ks
        self.gam_ks = gam_ks
        if(0):
            print(f'LowTcoef: {self.lowTcoef}')
            print(f'C_ks: {C_ks}')
            print(f'gam_ks: {gam_ks}')
            sys.exit(0) #exit the program after printing the coefficients
        return 

    # output TCF for a given set of C_ks and gam_ks
    def TCF(self,plotme=False,ax=plt,mode=None):
        if mode is None: mode = self.mode #allowing override of the mode from the __init__

        self.get_coeffs(mode=mode)

        t = np.linspace(0,5,1000)
        C = np.zeros_like(t,dtype=complex)
        # C += self.C_ks[0]*np.exp(-self.gam_ks[0]*t)
        print(f'C0: {self.C_ks[0]}, mode : {mode}')
        for k in range(1,self.N_exp):
            C += self.C_ks[k]*np.exp(-self.gam_ks[k]*t)
        if plotme:
            ax.plot(t,C.real)
            print(mode)
            plt.show()
        return t,C

    # Calculate the C_ks and gam_ks for a bath, mode is the bath decomposition mode
    def get_coeffs(self,mode=None):
        if mode is None: mode = self.mode #allowing override of the mode from the __init__

        if mode == 'highT':
            return np.array([self.C0hot]),np.array([self.gam]) #the high temperature limit - no matsubara terms

        if mode in ['matsubara','nmats']: # Generate the K matsubara frequencies (the nmats will have the same freqs, but different c0 and no low temp truncation)
            betaN = self.beta/self.N_mds
            wN=1/(betaN*self.hbar)
            wns = np.array([2*wN*np.pi*k/self.N_mds for k in range(0,self.mu+1)])
            print(f'mode {mode} with {self.mu} pairs of matsubara modes')
            # print(wns)
            return self.calc_coefs(wns)

        if mode == 'nbead':
            betaN =  self.beta/self.N_mds
            wN=1/(betaN*self.hbar)
            wks = np.array([2*wN*np.sin(np.pi*k/self.N_mds) for k in range(0,self.mu+1)])
            print(f'mode {mode} with {self.mu} pairs of matsubara modes')
            # print(wks)
            return self.calc_coefs(wks)
        else:
            raise ValueError('Invalid mode')
        return

    def get_C_UDs(self):
        # Calculate the coefficients C_U, c_D_LEFT, c_D_RIGHT for the bath (that are used in the FEOM code)
        self.c_U = np.zeros((self.N_exp,self.L+1),dtype=complex)
        self.c_D_LEFT = np.zeros((self.N_exp,self.L+1),dtype=complex)
        self.c_D_RIGHT = np.zeros((self.N_exp,self.L+1),dtype=complex)
        for ki in range(self.N_exp):
            for nk in range(self.L+1):
                self.c_U[ki,nk] = np.sqrt((nk+1)*abs(self.C_ks[ki]))
                if abs(self.C_ks[ki]) < 1e-10:
                    self.c_D_LEFT[ki,nk] = 0.0
                    self.c_D_RIGHT[ki,nk] = 0.0
                    print(f'Warning: C_ks({ki}) is zero, setting superoperator terms to zero')
                else:
                    self.c_D_LEFT[ki,nk] = -np.sqrt(nk/abs(self.C_ks[ki]))*self.C_ks[ki]
                    self.c_D_RIGHT[ki,nk] = np.sqrt(nk/abs(self.C_ks[ki]))*np.conj(self.C_ks[ki])
        return 



