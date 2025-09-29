import numpy as np
import matplotlib.pyplot as plt
import sys,copy
import Feom.baths.utils as utils

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
        self.LTCorr = getattr(params, 'LTCorr', None)      
        self.bathmode = params.bathmode
        self.L = params.L                      # max tier of the ADOs
        self.K = params.K                      # number of exponential terms in the BCF
        # General parameters
        self.eta = params.eta
        self.gam = params.gam
        self.beta = params.beta
        self.hbar = params.hbar
        # Paramaters for bath indexing
        self.N_nonmats = 1                          # number of exponential modes in BCF that are NOT matsubara terms (the Temp. ind. exp.)
        self.mu = params.K                          # number of pairs of matsubara modes/ r.p. modes, each pair gives a single exponential term
        self.N_exp = self.N_nonmats + self.mu       # number of exponential terms in the BCF [Temp ind. Exponential, <--- Matsubara Exponentials --->]
        self.N_exp_prop = self.N_nonmats + self.mu  # number of exponentials that we are propogating (should be the same as N_exp unless we want to truncate some matsubara terms (IT not included here))
        self.N_mds = 2*self.mu+1                    # number of individual beads or matsubara modes (ODD)
        #
        self.mode= params.bathmode
        ### Calculate the C_ks and gam_ks for the bath and add to the class
        self.get_coeffs()

    def J(self,w,plotme=False,ax=plt):
        w = np.linspace(0,2,1000)
        Jw = self.eta*self.gam*w/(w**2+self.gam**2)
        if plotme : 
            ax.plot(w,Jw)
        return Jw

    # Calculate the C_ks and gam_ks for a given set of ws
    def calc_coefs(self):
        d = np.zeros(self.N_exp)
        ### Calculate the 0th (non-matsubara) term
        if self.mode in ['nbead','nmats']:              # give the finite mode prefactor
            d0sum = np.sum(1/(self.gam**2*np.ones_like(self.ws[1:]) - self.ws[1:]**2))
            d[0] = (self.hbar*self.eta*self.gam/2) * (2/(self.beta*self.hbar*self.gam) + (4*self.gam/(self.beta*self.hbar))*d0sum)
        elif self.mode == 'matsubara': # give the infinite mode prefactor
            if self.LTCorr in ['iIT','viIT']:    # improved Ishizaki-Tanimura termination changes d0 (1/wn^2 instead of 1/(wn^2-gam^2))
                d0sum = (self.gam**2)*np.sum(1/((self.gam**2*np.ones_like(self.ws[1:]) - self.ws[1:]**2)*self.ws[1:]**2))
                d[0] =  (self.hbar*self.eta*self.gam/2) * (2/(self.beta*self.hbar*self.gam) + (4*self.gam/(self.beta*self.hbar))*d0sum-self.gam*self.beta*self.hbar/6)
            else: # standard ishizaki-tanimura
                d[0] = (self.hbar*self.eta*self.gam/2) /np.tan(self.beta*self.hbar*self.gam/2)
        else:
            raise ValueError('Invalid mode')
        ### Calculate the rest of the terms (matsubara terms)
        d[1:] = -(2*self.eta*self.gam/self.beta) * self.ws[1:]/(self.gam**2*np.ones_like(self.ws[1:]) - self.ws[1:]**2)
        dI = -self.hbar*self.eta*self.gam/2

        # Calculate the C_ks
        self.C_ks = np.zeros(self.N_exp,dtype=complex)
        self.C_ks[0] = (d[0] + dI*1.j)
        self.C_ks[1:] = d[1:]

        #Calculate the gam_ks
        self.gam_ks = np.zeros(self.N_exp,dtype=complex)
        self.gam_ks[0] = self.gam
        self.gam_ks[1:] = self.ws[1:]
        
        # Calculate the low temperature coefficient  (for use in baths/utils.py)
        if self.mode == 'matsubara':  # The Ishizaki-Tanimura terminator coefficient
            if self.LTCorr in ['iIT','viIT']:
                fac = -2*self.gam*self.eta/(self.beta*self.hbar**2)
                summ = np.sum(1/(self.ws[1:]**2))
                self.lowTcoef =  fac*(self.beta**2*self.hbar**2/24 - summ)
            else:
                self.lowTcoef = self.eta* ((1/(2*self.hbar))* ((1/(np.tan(self.beta*self.hbar*self.gam/2))) - (2/(self.beta*self.hbar*self.gam))))  ### Terms without removing of the matsubara terms that have been included
                self.lowTcoef = self.lowTcoef -  self.eta*(2*self.gam/(self.beta*self.hbar**2))*np.sum(1/(self.gam**2*np.ones_like(self.ws[1:]) - self.ws[1:]**2))  ### remove the Matsubara terms that have been explicitly included
        
        # Calculate constant term k for viIT
        if self.mode == 'matsubara' and self.LTCorr == 'viIT':
            ''' Calculae the variationally fitted k term for the BCF.
            This uses the fitting of Rg function, so needs cothpoles backend.
            Currently done in a hacky way by setting up new bath object with the same params'''
            from Feom.baths.debye.debye_cothpoles import Debye_cothpoles
            params=copy.deepcopy(self)
            params.bathmode='Pade[N/N]'
            params.LTCorr = None
            params.save_debug_data = False
            params.plot_debug_data = False
            padebath = Debye_cothpoles(params)
            Pw,support  = padebath.get_support_and_values()
            
            ### Calculate the Rg corresponding to the current set of matsubara frequencies
            wns = self.ws[1:]
            Pw_iIT=np.zeros_like(support) # Pole function for the iIT
            for wn in wns:
                Pw_iIT += 1/(support**2+wn**2)
                Pw_iIT -= 1/(wn**2)  # remove from the infinite sum 1/wn^2 (below)
            Pw_iIT += (self.beta*self.hbar)**2/24 #the infinite sum correction

            # Now calculate the k term from the difference in the two pole functions
            diff = Pw - Pw_iIT
            # now find the best fit constant to this difference
            self.k = np.mean(diff)
            if(0):
                plt.plot(support,Pw,label='exact')
                plt.plot(support,Pw_iIT,label='iIT')
                plt.plot(support,diff,label='diff')
                plt.plot(support,Pw_iIT+self.k,label='adjusted, k={:.2e}'.format(self.k))
                plt.legend()
                plt.show()
            np.savetxt('iIT_k.txt',[self.k])

            # add on the contribution to d0
            self.C_ks[0]+= - 2*self.eta*self.gam**2*self.k/self.beta        

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

        # if mode == 'highT':
        #     return np.array([self.C0hot]),np.array([self.gam]) #the high temperature limit - no matsubara terms

        if mode in ['matsubara','nmats']: # Generate the K matsubara frequencies (the nmats will have the same freqs, but different c0 and no low temp truncation)
            betaN = self.beta/self.N_mds
            wN=1/(betaN*self.hbar)
            wns = np.array([2*wN*np.pi*k/self.N_mds for k in range(0,self.mu+1)])
            # print(f'mode {mode} with {self.mu} pairs of matsubara modes')
            # print(wns)
            self.ws = wns
            return self.calc_coefs()

        if mode == 'nbead':
            betaN =  self.beta/self.N_mds
            wN=1/(betaN*self.hbar)
            wks = np.array([2*wN*np.sin(np.pi*k/self.N_mds) for k in range(0,self.mu+1)])
            # print(f'mode {mode} with {self.mu} pairs of matsubara modes')
            # print(wks)
            self.ws = wks
            return self.calc_coefs()
        else:
            raise ValueError('Invalid mode')
        return

   



