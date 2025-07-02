import numpy as np
import matplotlib.pyplot as plt
import sys

# A class for the Debye bath
''' A class to represent the Debye bath for the FEOM code, with AAA decomposition of A(w)/W

eta : Coupling strength
gam : Cuttoff frequency

The spectral density is given by:
J(w) = \frac{\eta\gamma\omega}{\omega^2 + \gamma^2}
using the discretized definition, giving coefficients:
J(w) = (\pi/2) * \sum_{\alpha} \frac{c_\alpha^2}{m_\alpha \omega_\alpha} \delta(w - w_\alpha)
'''

# We need to add in pade approximants, but otherwiseshould be mostly complete

class Debye_colepoles():
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
        self.N_nonmats = 1                      # number of exponential modes in BCF that are NOT matsubara terms (the Temp. ind. exp.)
        self.mu = params.K                      # number of pairs of matsubara modes/ r.p. modes, each pair gives a single exponential term
        self.N_exp = self.N_nonmats + self.mu   # number of exponential terms in the BCF [Temp ind. Exponential, <--- Matsubara Exponentials --->]
        self.N_mds = 2*self.mu+1                # number of individual beads or matsubara modes (ODD)
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
    
    def P(self,w): # Pole function for the coth, that we are gonna approximate
        return (1/w)*(1/np.tanh(self.beta*self.hbar*w/2) - 2/(self.beta*self.hbar*w))  
    
    def calc_coefs(self):
        support = np.linspace(-1000,1000,10000,dtype=np.complex128) # support for the AAA decomposition
        values = self.S(support)                 # values of the spectral density at the support points
        # Use the AAA decomposition to get the coefficients
        # save the support and values to be read by matlab
        folder = f'{script_dir}/aaa/.files'  # folder to save the data to
        data = np.column_stack((support.real, values.real, values.imag))  
        np.savetxt(f'{folder}/aaa_data.txt', data, header='Support Re[Values] Im[Values]', comments='')
        print('Running AAA decomposition in MATLAB...')
        # run using os
        # os.system(f"matlab -batch 'cd {script_dir}/aaa;run_aaa_fromfile({self.N_exp//2})' > /dev/null 2>&1") # run the matlab script to get the AAA coefficients
        print('AAA decomposition complete, loading results...')
        ### Load the aaa results
        repoles = np.loadtxt(f'{folder}/pol_real.txt')
        impoles = np.loadtxt(f'{folder}/pol_imag.txt')
        self.poles = repoles + 1.j * impoles
        reres = np.loadtxt(f'{folder}/res_real.txt')
        imres = np.loadtxt(f'{folder}/res_imag.txt')
        self.res = reres + 1.j * imres
        if self.N_exp != len(self.poles):
            print(f'Number of poles {len(self.poles)} does not match number of residues {len(self.res)}, changing now')
        self.N_exp = len(self.poles)  # update the number of exponentials

    # Calculate the C_ks and gam_ks for a given set of ws
    def calc_coefs_old(self,ws):
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

        # Calculate the low temperature coefficient for LowT correction
        # self.lowTcoef = self.eta/(self.beta*self.hbar**2) - (1/self.hbar**2)*np.sum(np.real(C_ks)/gam_ks) if self.bathmode == 'matsubara' else 0 #OLD ONE
        
        # Calculate the low temperature coefficient for LowT correction ### NEW ONE
        self.lowTcoef = self.eta* ((1/(2*self.hbar))* ((1/(np.tan(self.beta*self.hbar*self.gam/2))) - (2/(self.beta*self.hbar*self.gam))))  ### Terms without removing of the matsubara terms that have been included
        self.lowTcoef = self.lowTcoef -  self.eta*(2*self.gam/(self.beta*self.hbar**2))*np.sum(1/(self.gam**2*np.ones_like(ws[1:]) - ws[1:]**2))  if self.mode == 'matsubara' else 0 ### remove the Matsubara terms that have been explicitly included

        self.C_ks = C_ks
        self.gam_ks = gam_ks
        if(0):
            print(f'LowTcoef: {self.lowTcoef}')
            print(f'C_ks: {C_ks}')
            print(f'gam_ks: {gam_ks}')
            sys.exit(0) #exit the program after printing the coefficients
        return 
        # return C_ks,gam_ks

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



