import numpy as np
import matplotlib.pyplot as plt
import sys, os
import Feom.baths.coth_decomp.cothPade as cothPade 
import Feom.baths.coth_decomp.cothAAA as cothAAA
import Feom.baths.utils as utils

# A class for the Debye bath
''' A class to represent the Debye bath for the FEOM code, with AAA/Pade decomposition of A(w)/w

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
        self.terminate = hasattr(params,'LTCorr') # whether or not we are using a terminator
        # General parameters
        self.eta = params.eta
        self.gam = params.gam
        self.beta = params.beta
        self.hbar = params.hbar
        # Paramaters for bath indexing
        self.N_nonmats = 1                      # number of exponential modes in BCF that are NOT matsubara terms (the Temp. ind. exp.)
        self.mu = params.K                      # number of pairs of matsubara modes/ r.p. modes, each pair gives a single exponential term
        self.N_exp_prop = self.N_nonmats + self.mu   # number of exponential terms in the BCF EXPLICITLY PROPOGATED [Temp ind. Exponential, <--- Matsubara Exponentials --->]
        #
        self.mode= params.bathmode
        # self.C0hot = self.eta/self.beta -1.j*self.hbar*self.eta*self.gam/2 # C_0 with no matsubara terms
        ### Calculate the C_ks and gam_ks for the bath and add to the class
        self.get_coefs()
        # self.TCF(plotme=True,ax=plt) # Calculate the TCF for the bath and plot it
        # params.K = len(self.w_i) # Update the number of exponentials in the params object to match the number of poles found
        # params.lowTCorr=True if self.k != 0 else False


    def J(self,w,plotme=False,ax=plt):
        Jw = self.eta*self.gam*w/(w**2+self.gam**2)
        if plotme : 
            ax.plot(w,Jw)
        return Jw
    
    def P(self,w): # Pole function for the coth, that we are gonna approximate
        return (self.beta*self.hbar/4)*(1/w)*(1/np.tanh(self.beta*self.hbar*w/2) - 2/(self.beta*self.hbar*w))  
    
    def P_aaa(self,w): # Pole function for the coth, calculated using the original poles and residues from the AAA decomposition
        result = 0.0 if np.isscalar(w) else np.zeros_like(w,dtype=np.complex128)
        for k in range(len(self.res_original)):
            result += self.res_original[k] / (w - self.poles_original[k])
        return result
    
    def P_aaa_realcoeffs(self,w): # Pole function for the coth, calculated using the poles and residues from the AAA decomposition, with real coefficients
        result = 0.0 if np.isscalar(w) else np.zeros_like(w,dtype=np.complex128)
        for k in range(len(self.gam_i)):
            result += self.gam_i[k] / (w**2 + self.w_i[k]**2)
        return result
         
    
    def calc_poles(self): # Reclusters the poles and residues from coth(x)
        '''
        self.terminate = TRUE/FALSE is the switch for whether or not we are using a terminator
        if there is no terminator, N_exp=N_exp_prop
        if not, N_exp=N_exp_prop + (enough terms to converge the low temp correction)
        the C_ks and gam_ks calculated will then be used in the terminator module to calcuate the terminator
        '''

        if self.bathmode=='AAA':
            print(f'Using AAA decomposition for the bath.')
            if(0): # Calculate the support of the coth function such that J(w)/w is sampled evenly [DOESNT WORK WELL]
                N_support = 100000
                x_j=np.arange(1,N_support+1)/(self.gam**2*(N_support+1)) #equally spaced x
                support = np.concatenate([np.array([-200]),
                    -np.abs(np.sqrt(1 / x_j - self.gam**2)),
                    np.flip(np.abs(np.sqrt(1 / x_j - self.gam**2)))
                    ,np.array([200])])
                values= self.P(support) # values of the coth function at the support points
                self.gam_i, self.w_i, self.k, mu_tot = cothAAA.get_coeffs(self,support,values,self.terminate)
            else: # Use support with uniform spacing in w (done within the cothAAA module)
                self.gam_i, self.w_i, self.k, mu_tot = cothAAA.get_coeffs(self,None,None,self.terminate) 

        elif self.bathmode[0:4]=='Pade':
            Padetype = self.bathmode[4:] # get the type of Pade decomposition
            print(f'Using Pade decomposition of type {Padetype} for the bath.')
            eta, xi, R_N, mu_tot =cothPade.get_coeffs(Padetype,self.mu,self.terminate) # get the poles and residues from the pade module
            ### convert to the same format as the AAA decomposition
            self.gam_i = eta
            self.w_i = xi/(self.beta *self.hbar)
            self.k = R_N*(self.beta*self.hbar)**2/2.
        else:
            raise ValueError('Invalid tyre of coth decomposition specified. Use "AAA" or "Pade..." .')
        
        self.N_exp = self.N_nonmats + mu_tot  # total number of exponentials in the BCF 

        if(1): # save and Plot the approximated function and J(w)
            w = np.linspace(-250,250,20000,dtype=np.complex128) 
            values = self.P(w)                     # values of the pole function at the w points
            plt.figure(figsize=(5,5))
            plt.plot(w.real, values.real, label='Original Function', color='blue')
            plt.plot(w.real, self.P_aaa_realcoeffs(w).real+self.k, label=f'{self.bathmode} Approximation, imaginary poles/residues', color='green')
            # plt.plot(w, self.P_aaa(w).real, label='AAA Approximation', color='red')
            plt.plot(w.real, self.J(w).real, label='J(w)', color='orange')
            plt.xlabel('w')
            plt.ylabel('Function Value')
            plt.title(f'{self.mu} mode approximation of {self.bathmode} Approximation of the Pole Function')
            plt.legend()
            plt.grid()
            plt.show() if self.plot_debug_data else None
            filename = f'{self.cleanbathmode}_Cothapproximation.txt'.replace('[N/N]','NoN').replace('[N-1/N]','Nm1oN')
            data = np.column_stack((w.real, self.P_aaa_realcoeffs(w).real+self.k, values.real, self.J(w).real))
            np.savetxt(filename, data, header='# w Re[P_approx(w)] Re[P(w)] Re[J(w)]', comments='') if self.save_debug_data else None
            print(f'Saved the approximation plot to {filename}')
            # sys.exit(0) # exit the program after plotting the approximation
        if(1): #save and plot the poles and the Matsubara terms if they were to be used
            wmax= np.max(np.abs(self.w_i)) 

            fig, pole_ax = plt.subplots(figsize=(5,5))
            pole_ax.set_title(r'Poles and residues: $\sum_i \frac{\gamma_i}{\omega^2+\omega_i^2}$', fontsize=16)
            pole_ax.set_xlabel(r'$\gamma_i$', fontsize=14)
            pole_ax.set_ylabel(r'$\omega_i$', fontsize=14, rotation=0, labelpad=10)
            pole_ax.grid(True)
            pole_ax.plot(self.gam_i, self.w_i, 'x', label=f'AAA, N={len(self.w_i)}', color='blue')
            pole_ax.legend()
            
            data = np.column_stack((self.gam_i.real, self.w_i.real))
            filename= f'{self.cleanbathmode}_poles_and_residues.txt'.replace('[N/N]','NoN').replace('[N-1/N]','Nm1oN')
            np.savetxt(filename, data, header=f'# gamma_i w_i N={len(self.w_i)}', comments='') if self.save_debug_data else None
            if self.bathmode in ['Pade[N/N]','AAA']:
                k_filename= f'{self.cleanbathmode}_k.txt'.replace('[N/N]','NoN').replace('[N-1/N]','Nm1oN')
                np.savetxt(k_filename, [self.k], header=f'# k', comments='') if self.save_debug_data else None

            plt.show() if self.plot_debug_data else None


    # Calculate the C_ks and gam_ks
    def get_coefs(self):
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

        #Calculate the gam_ks 
        self.gam_ks = np.zeros(self.N_exp,dtype=complex)
        self.gam_ks[0] = self.gam
        self.gam_ks[1:] = self.w_i[:]

        # Calculate the low temperature coefficient
        self.lowTcoef=-2*self.eta*self.gam*self.k/(self.beta*self.hbar**2) 

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
            ax.plot(t,C.real,label='AAA tcf real')
            ax.plot(t,C.imag,label='AAA tcf imag')
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
    





