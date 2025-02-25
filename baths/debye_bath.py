import numpy as np
import matplotlib.pyplot as plt

# A class for the Debye bath

# We need to add in pade approximants, but otherwiseshould be mostly complete

class Debye_bath():
    def __init__(self,params):
        # General parameters
        self.eta = params.eta
        self.gam = params.gam
        self.beta = params.beta
        self.hbar = params.hbar
        # Paramaters for bath indexing
        self.N_nonmats = 1                      # number of exponential modes in BCF that are NOT matsubara terms (the Temp. ind. exp.)
        self.mu = params.K                      # number of pairs of matsubara modes/ r.p. modes, each pair gives a single exponential term
        self.N_exp = self.N_nonmats + self.mu   # number of exponential terms in the BCF [Temp ind. Exponential, <--- Matsubara Exponentials --->]
        self.N_mds = 2*self.mu+1                # number of individual beads or matsubara modes (ODD)
        #
        self.mode= params.bathmode
        self.C0hot = self.eta/self.beta -1.j*self.hbar*self.eta*self.gam/2 # C_0 with no matsubara terms

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
        if self.mode == 'nbead':
            d0sum = np.sum(1/(self.gam**2 - ws**2))
            d[0] = self.eta/self.beta + (2*self.eta*self.gam**2/self.beta) * d0sum 
        elif self.mode == 'matsubara': #give the infinite mode prefactor
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


        return C_ks,gam_ks

    # output TCF for a given set of C_ks and gam_ks
    def TCF(self,plotme=False,ax=plt,mode=None):
        if mode is None: mode = self.mode #allowing override of the mode from the __init__

        C_ks,gam_ks = self.get_coeffs(mode=mode)

        t = np.linspace(0,100,1000)
        C = np.zeros_like(t,dtype=complex)
        C += C_ks[0]*np.exp(-gam_ks[0]*t)
        for k in range(1,self.N_exp):
            C += C_ks[k]*np.exp(-gam_ks[k]*t)
        if plotme:
            ax.plot(t,C.real)
            print(mode)
        return t,C

    # Calculate the C_ks and gam_ks for a bath, mode is the bath decomposition mode
    def get_coeffs(self,mode=None):
        if mode is None: mode = self.mode #allowing override of the mode from the __init__

        if mode == 'highT':
            return np.array([self.C0hot]),np.array([self.gam]) #the high temperature limit - no matsubara terms

        if mode == 'matsubara':
            betaN = self.beta/self.N_mds
            wN=1/(betaN*self.hbar)
            wns = np.array([2*wN*np.pi*k/self.N_mds for k in range(0,self.mu+1)])
            # print(wns)
            return self.calc_coefs(wns)

        if mode == 'nbead':
            betaN =  self.beta/self.N_mds
            wN=1/(betaN*self.hbar)
            wks = np.array([2*wN*np.sin(np.pi*k/self.N_mds) for k in range(0,self.mu+1)])
            # print(wks)

            return self.calc_coefs(wks)
        else:
            raise ValueError('Invalid mode')
        return





