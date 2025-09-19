import numpy as np
import matplotlib.pyplot as plt
import sys, os
import Feom.baths.coth_decomp.cothPade as cothPade 
import Feom.baths.coth_decomp.cothAAA as cothAAA
import Feom.baths.utils as utils

# A class for the Debye bath
''' A class to represent the Debye bath for the FEOM code, 
with bathmode options:
Pade... : Pade decomposition of A(w)/w
AAA     : AAA decomposition of A(w)/w
AAAmc   : AAA decomposition of A(w)/w with high frequency modes treated markovianly

eta : Coupling strength
gam : Cuttoff frequency

The spectral density is given by:
J(w) = \frac{\eta\gamma\omega}{\omega^2 + \gamma^2}
using the discretized definition, giving coefficients:
J(w) = (\pi/2) * \sum_{\alpha} \frac{c_\alpha^2}{m_\alpha \omega_\alpha} \delta(w - w_\alpha)
'''


class Debye_cothpoles():
    def __init__(self,params):
        self.save_debug_data = True
        self.plot_debug_data = False
        self.bathmode = params.bathmode
        self.cleanbathmode = self.bathmode.replace(' ','_').replace('/','_') # clean the bath mode name for saving files
        self.L = params.L                      # max tier of the ADOs
        #NOTE the above was used under the assumption that termination was of the extra terms in the BCFs (for the FAY way self.N_exp_prop=self.N_exp always)
        # NEED TO LOOK INTO THIS ^^ maybe this is the best way, but IDK ^^
        # General parameters
        self.eta = params.eta
        self.gam = params.gam
        self.beta = params.beta
        self.hbar = params.hbar
        # Paramaters for bath indexing
        self.N_nonmats = 1                      # number of exponential modes in BCF that are NOT matsubara terms (the Temp. ind. exp.)
        self.mu = params.K                      # number of pairs of matsubara modes/ r.p. modes, each pair gives a single exponential term
        #
        self.mode= params.bathmode
        ### Calculate the C_ks and gam_ks for the bath and add to the class
        self.get_coefs(params)
        print(f'Calculated {self.N_exp} exponentials for the bath with {self.mu} coth poles')
        # self.TCF(plotme=True,ax=plt) # Calculate the TCF for the bath and plot it


    def J(self,w,plotme=False,ax=plt):
        Jw = self.eta*self.gam*w/(w**2+self.gam**2)
        if plotme : 
            ax.plot(w,Jw)
        return Jw
    
    def P(self,w,M=1): 
        ''' P(w), the pole function for the coth.
        M is the number of low frequency matsubara terms to approximate by the highest frequency term'''
        Pw=np.zeros_like(w)
        for i, wi in enumerate(w):
            if wi!=0:
                Pw[i] = (self.beta*self.hbar/4)*(1/wi)*(1/np.tanh(self.beta*self.hbar*wi/2) - 2/(self.beta*self.hbar*wi))
            else: # treat the 0 divergence nicely
                Pw[i] = ((self.beta*self.hbar)**2)/24
        # Add on the M approximated low frequency terms, if needed
        wn_s = np.array([2*np.pi*n/(self.beta*self.hbar) for n in range(1,M+1)])
        for wn in wn_s:
            Pw -= 1/(w**2 + wn**2)
        Pw += M/(w**2+wn_s[-1]**2) # add back on the terms (with all of the same frequencies)
        return Pw
    
    def P_aaa(self,w):
        '''Calculate P(w) using the COMPLEX poles from the AAA decomposition of coth(x) [no k added]'''
        result = 0.0 if np.isscalar(w) else np.zeros_like(w,dtype=np.complex128)
        for k in range(len(self.res_original)):
            result += self.res_original[k] / (w - self.poles_original[k])
        return result
    
    def P_aaa_realcoeffs(self,w): 
        ''' 'Calculate P(w) using the IMAGINARY poles from the AAA decomposition of coth(x) [no k added]'''
        result = 0.0 if np.isscalar(w) else np.zeros_like(w,dtype=np.complex128)
        for k in range(len(self.gam_i)):
            result += self.gam_i[k] / (w**2 + self.w_i[k]**2)
        return result
         
    def get_support_and_values(self, mode='uniform',N_support = 100000):
        ''' Generate the support and values for the AAA decomposition of the pole function'''    
        if mode=='log': # logarithmic spacing including zero
            eps=1e-4
            w_max=self.gam # start with the cuttoff frequency
            Jw_min_tol = 1e-4      # tolerance for the maximum frequency of the grid for the AAA decomposition
            while self.J(w_max) > Jw_min_tol: w_max += 10 # find the maximum frequency where J(w) is still non-zero
            x_pos = np.logspace(np.log10(eps), np.log10(w_max), N_support // 2)
            support = np.concatenate((-x_pos[::-1], [0.0], x_pos))
            self.support_param_str = f'log_N{N_support}_wmax{int(w_max)}'
        elif mode=='quadrature': # use the points from the quadrature of J(w)
            x_j = np.arange(1, N_support + 1)/((N_support + 1)*(self.gam**2))  
            w_j = np.sqrt(1/x_j - self.gam**2)  # abscissas
            support = np.concatenate((-w_j[::-1], [0.0], w_j))  # support points
            self.support_param_str = f'gauss_legendre_N{N_support}_gam{self.gam}' # save the parameters used to generate the support points
        elif mode == 'uniform': # uniform spacing including zero
            w_max=self.gam # start with the cuttoff frequency
            Jw_min_tol = 1e-5      # tolerance for the maximum frequency of the grid for the AAA decomposition
            while self.J(w_max) > Jw_min_tol: w_max += 10
            w_max=200
            support = np.linspace(-w_max,w_max,N_support,dtype=np.complex128)
            self.support_param_str = f'uniform_N{N_support}_wmax{int(w_max)}' # save the parameters used to generate the support points
        elif mode == 'arctanh': #NOTE - need to choose w_max carefully here
            w_max=200
            eps=1e-5
            range = np.linspace(-1+eps, 1-eps, N_support)
            x = np.arctanh(range) * w_max
            support = x[1:-1]
            print(support)
            self.support_param_str = f'arctanh_N{N_support}_wmax{int(w_max)}'
        else:
            raise ValueError('Invalid mode for generating support points. Use "log", "quadrature" or "uniform".')

        values= self.P(support) # values of the POLE function at the support points
        return values, support
    
    def calc_poles(self): 
        '''recluster poles from the coth function either using Pade or AAA algos'''
        if self.bathmode[0:3]=='AAA':
            print(f'Using AAA decomposition for the bath.')
            values, support = self.get_support_and_values()
            Rg = values*(2/self.beta)   # so that we fit the Rg (not the pole function directly)
            Rg_gam_i, self.w_i, Rg_k, mu_tot = cothAAA.get_coeffs(self,support,Rg)
            self.gam_i = Rg_gam_i*(self.beta/2) # convert back to gam_i (as we originally fitted pole function in earlier code)
            self.k = Rg_k*(self.beta/2)         # convert back to k

        elif self.bathmode[0:4]=='Pade':
            Padetype = self.bathmode[4:] # get the type of Pade decomposition
            print(f'Using Pade decomposition of type {Padetype} for the bath.')
            eta, xi, R_N, mu_tot =cothPade.get_coeffs(Padetype,self.mu) # get the poles and residues from the pade module
            ### convert to the same format as the AAA decomposition
            self.gam_i = eta
            self.w_i = xi/(self.beta *self.hbar)
            self.k = R_N*(self.beta*self.hbar)**2/2.
        else:
            raise ValueError('Invalid type of coth decomposition specified. Use "AAA" or "Pade..." .')
        
        self.N_exp = self.N_nonmats + mu_tot  # total number of exponentials in the BCF 


    # Calculate the C_ks and gam_ks
    def get_coefs(self,params):
        self.calc_poles()  # Calculate the poles and residues from the coth function, and set the w_i and gam_i attributes

        d = np.zeros(self.N_exp)
        ### Calculate the 0th term
        d0sum = np.sum(self.gam_i/(self.gam**2 - self.w_i**2))
        d[0] = self.eta/self.beta + (2*self.eta*self.gam**2/self.beta)*d0sum - 2*self.eta*self.gam**2*self.k/self.beta
        # NOTE: TOM FAYS CODE DOES NOT INCLUDE k IN C0. FOR THIS REASON THE PADE[N/N] DOES NOT AGREE WITH HIS RESULTS
        # IF k IS REMOVED, THEN THE RESULTS AGREE.

        ### Calculate the rest of the terms (matsubara terms)
        d[1:] = -(2*self.eta*self.gam/self.beta) * self.w_i[:]*self.gam_i[:]/(self.gam**2*np.ones_like(self.w_i) - self.w_i[:]**2)
        dI = -self.hbar*self.eta*self.gam/2

        # Calculate the C_ks
        self.C_ks = np.zeros(self.N_exp,dtype=complex)
        self.C_ks[0] = (d[0] + dI*1.j)
        self.C_ks[1:] = d[1:]

        # Calculate the gam_ks 
        self.gam_ks = np.zeros(self.N_exp,dtype=complex)
        self.gam_ks[0] = self.gam
        self.gam_ks[1:] = self.w_i[:]


        # Remove high frequencies if using AAAmc
        if self.bathmode=='AAAmc':
            cothAAA.markovian_pole_trunc(self)
            # Update K in parameters (changed for the bath in markovian_pole_trunc)
            params.K = self.K

        # Printouts and debug data
        w = np.logspace(1e-10, 100, 10000,dtype=np.complex128)  # 200000 points from 1e-10 to 100
        w = np.concatenate((-w[::-1], [0.0], w))  # support points with 0
        w = np.concatenate([np.linspace(0,2,20000,dtype=np.complex128) ,np.linspace(2,200,500,dtype=np.complex128)])
        values = self.P(w)  # values of the pole function at the w points          
        if(self.plot_debug_data): 
            ### Plot the approximated function and J(w)
            fit_fig,fit_ax=plt.subplots(figsize=(5,5))
            fit_ax.plot(w.real, values.real, label='Original Function', color='blue')
            fit_ax.plot(w.real, self.P_aaa_realcoeffs(w).real+self.k, label=f'{self.bathmode} Approximation, imaginary poles/residues', color='green')
            fit_ax.plot(w.real, self.J(w).real, label='J(w)', color='orange')
            fit_ax.set_xlabel('w')
            fit_ax.set_ylabel('Function Value')
            fit_ax.set_title(f'{self.mu} mode approximation of {self.bathmode} Approximation of the Pole Function')
            fit_ax.legend()
          
            #save and plot the poles and the Matsubara terms if they were to be used
            fig, pole_ax = plt.subplots(figsize=(5,5))
            pole_ax.set_title(r'Poles and residues: $\sum_i \frac{\gamma_i}{\omega^2+\omega_i^2}$', fontsize=16)
            pole_ax.set_xlabel(r'$\gamma_i$', fontsize=14)
            pole_ax.set_ylabel(r'$\omega_i$', fontsize=14, rotation=0, labelpad=10)
            pole_ax.grid(True)
            pole_ax.plot(self.gam_i, self.w_i, 'x', label=f'AAA, N={len(self.w_i)}', color='blue')
            pole_ax.legend()
            plt.show()
        if(self.save_debug_data):
            # save the approximatoin of the pole function
            filename = f'{self.cleanbathmode}_Cothapproximation.txt'.replace('[N/N]','NoN').replace('[N-1/N]','Nm1oN')
            data = np.column_stack((w.real, self.P_aaa_realcoeffs(w).real+self.k, values.real, self.J(w).real))
            np.savetxt(filename, data, header='# w Re[P_approx(w)] Re[P(w)] Re[J(w)]', comments='')
            print(f'Saved the approximation plot to {filename}')
            # save the poles and residues
            data = np.column_stack((self.gam_i.real, self.w_i.real))
            filename= f'{self.cleanbathmode}_poles_and_residues.txt'.replace('[N/N]','NoN').replace('[N-1/N]','Nm1oN')
            np.savetxt(filename, data, header=f'# gamma_i w_i N={len(self.w_i)}', comments='') 
            if self.bathmode in ['Pade[N/N]','AAA']:
                k_filename= f'{self.cleanbathmode}_k.txt'.replace('[N/N]','NoN').replace('[N-1/N]','Nm1oN')
                np.savetxt(k_filename, [self.k], header=f'# k', comments='') 

        return 

    # output TCF for a given set of C_ks and gam_ks
    def TCF(self,plotme=False,ax=plt,mode=None):
        if mode is None: mode = self.mode #allowing override of the mode from the __init__
        def C_analytic(t):
            C_0 = (self.hbar*self.eta*self.gam/2)*(1/np.tan(self.beta*self.hbar*self.gam/2) - 1.j) 
            mu =1000
            M = 2*mu+1
            betaN = self.beta/M
            wN=1/(betaN*self.hbar)
            wns = np.array([2*wN*np.pi*k/M for k in range(0,mu+1)])
            C_n = -(2*self.eta*self.gam/self.beta)* wns/(self.gam**2*np.ones_like(wns) - wns**2) 
            C_k= np.zeros_like(wns,dtype=complex)
            C_k[0] = C_0
            C_k[1:] = C_n[1:]
            gam_k = np.zeros_like(C_k,dtype=complex)
            gam_k[0] = self.gam
            gam_k[1:] = wns[1:]
            C = np.zeros_like(t,dtype=complex)
            for k in range(0,self.N_exp):
                C += C_k[k]*np.exp(-gam_k[k]*t)
            return C
        t = np.linspace(0,5,5000)
        C = np.zeros_like(t,dtype=complex)
        for k in range(0,len(self.C_ks)):
            C += self.C_ks[k]*np.exp(-self.gam_ks[k]*t)
        if plotme:
            ax.plot(t,C.real,label='tcf real')
            ax.plot(t,C.imag,label='tcf imag')
            C_analytic_t = C_analytic(t)
            ax.plot(t,C_analytic_t.real,'--',label='Analytic TCF real')
            ax.plot(t,C_analytic_t.imag,'--',label='Analytic TCF Imaginary')
            ax.xlabel('Time')
            ax.ylabel('TCF')
            ax.title(f'TCF for {self.bathmode} bath with {self.N_exp} exponentials')
            ax.legend()
            ax.grid(True)
            plt.show()
        return t,C
    





