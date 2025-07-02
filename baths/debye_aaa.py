import numpy as np
import matplotlib.pyplot as plt
import sys, os
import subprocess
# Get the directory of the current script file
script_dir = os.path.dirname(os.path.abspath(__file__))


# A class for the Debye bath
''' A class to represent the Debye bath for the FEOM code, with AAA decomposition.

NOTE: The parameter K is the total number of ADO indices, such that it MUST be even.
        This is so that there are equal numbers of l and m ADO indices.

eta : Coupling strength
gam : Cuttoff frequency

The spectral density is given by:
J(w) = \frac{\eta\gamma\omega}{\omega^2 + \gamma^2}
using the discretized definition, giving coefficients:
J(w) = (\pi/2) * \sum_{\alpha} \frac{c_\alpha^2}{m_\alpha \omega_\alpha} \delta(w - w_\alpha)
'''

# We need to add in pade approximants, but otherwiseshould be mostly complete

class Debye_aaa():
    def __init__(self,params):
        raise NotImplementedError('\nDebye_aaa is not implemented yet')
        # Bathmode and settings
        self.bathmode = params.bathmode
        self.L = params.L                      # max tier of the ADOs
        self.K = params.K                      # number of Matsubara terms in the BCF [Temp ind. Exponential, <--- Matsubara Exponentials --->]
        # General parameters
        self.eta = params.eta
        self.gam = params.gam
        self.beta = params.beta
        self.hbar = params.hbar
        # Paramaters for bath indexing
        self.N_nonmats = 0                      # set to 0 as the AAA will just give K/2 exponential terms in the BCF
        assert self.K % 2 == 0, "K must be even for AAA decomposition"
        self.N_exp = self.K//2                  # number of exponential terms in the BCF [Temp ind. Exponential, <--- Matsubara Exponentials --->]
        print(f'Using {self.N_exp} exponentials in the BCF for the bath {self.bathmode}')
        #
        self.mode= params.bathmode
        ### Calculate the C_ks and gam_ks for the bath and add to the class
        self.calc_coefs()
        self.TCF(plotme=True,ax=plt,mode=self.mode) # Calculate the TCF for the bath and plot it
        ### Calculate the coefficients C_U, c_D_LEFT, c_D_RIGHT for the bath (that are used in the FEOM code)
        self.get_C_UDs()

    def J(self,w,plotme=False,ax=plt):
        w = np.linspace(0,2,1000)
        Jw = self.eta*self.gam*w/(w**2+self.gam**2)
        if plotme : 
            ax.plot(w,Jw)
        return Jw,w

    def J_(self,w):
        return self.eta*self.gam*w/(w**2+self.gam**2)  # return the spectral density for a given w
    
    def A(self,w):
        return 1/(np.tanh(self.beta*self.hbar*w/2))  # return the A function for a given w [COULD PUT THE APPROX BITS IN HERE]
    
    def S(self,w): # return the spectral density for a given w
        fac = 2*np.pi*self.hbar
        return fac*self.J_(w)*(self.A(w)+1)  # return the spectral density for a given w        

    def S_aaa(self,w):
        """Return the spectral density for a given w using the AAA decomposition."""
        # Calculate the spectral density using the poles and residues
        S_aaa = np.zeros_like(w, dtype=complex)
        for k in range(len(self.poles)):
            S_aaa += self.res[k] / (w - self.poles[k])
        return S_aaa
    
    # Calculate the poles for the AAA decomposition
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


        # plot the poles
        if(0): # set to True to plot the poles
            plt.figure(figsize=(12, 6))
            pole_ax = plt.subplot(121)
            pole_ax.set_title('Poles')
            pole_ax.set_xlabel('Real Part')
            pole_ax.set_ylabel('Imaginary Part')
            pole_ax.grid(True)
            pole_ax.plot(repoles, impoles, 'o', label=f'AAA, N={len(repoles)}', color='blue')
            pole_ax.plot(reres, imres, 'x', label=f'Residues, N={len(reres)}', color='orange')
            pole_ax.legend()
            func_ax = plt.subplot(122)
            func_ax.set_title('Spectral Density')
            func_ax.set_xlabel('Frequency')     
            func_ax.set_ylabel('Spectral Density')
            func_ax.grid(True)
            w = np.linspace(-100, 100, 2000)  # Frequency range
            S_exact_values = np.array([self.S(omega) for omega in w])
            func_ax.plot(w, S_exact_values.real, label='Re[S_aaa(ω)]', color='red')
            func_ax.plot(w, S_exact_values.imag, label='Im[S_aaa(ω)]', color='green')
            S_aaa_values = self.S_aaa(w)
            func_ax.plot(w, S_aaa_values.real, label='Re[S_aaa(ω)]', color='purple', linestyle='--')
            func_ax.plot(w, S_aaa_values.imag, label='Im[S_aaa(ω)]', color='brown', linestyle='--')

            func_ax.legend()
            plt.tight_layout()
            plt.show()
        # define the support function for the AAA decomposition
        

    # output TCF for a given set of Residues and poles
    def TCF(self,plotme=False,ax=plt,mode=None):
        if mode is None: mode = self.mode #allowing override of the mode from the __init__

        fac = 1/(2*np.pi) 
        t = np.linspace(0,10,1000)
        C = np.zeros_like(t,dtype=complex)
        C_analytic = np.zeros_like(t,dtype=complex) # initialize the analytic TCF
        ws = np.linspace(-100,100,1000); dw=ws[1]-ws[0] # frequency range for the analytic TCF
        for i,ti in enumerate(t):
            C_analytic[i] = self.hbar*np.sum(self.J_(ws)*(self.A(ws)+1)*np.exp(-1.j*ws*ti))*dw # calculate the analytic TCF
        print('length of poles:',len(self.poles),'Nexp:',self.N_exp)
        C = np.zeros_like(t,dtype=complex) # initialize the TCF
        for k in range(0,self.N_exp):
            if self.poles[k].imag < 0:
                C += self.res[k]*np.exp(-1.j*self.poles[k]*t) 
        if plotme:
            ax.plot(t,C.real,label='AAA tcf real')
            ax.plot(t,C.imag,label='AAA tcf imag')
            ax.plot(t,C_analytic.real,'--',label='Analytic TCF real')
            ax.plot(t,C_analytic.imag,'--',label='Analytic TCF Imaginary')
            ax.xlabel('Time')
            ax.ylabel('TCF')
            ax.title(f'TCF for {self.bathmode} bath')
            ax.legend()
            ax.grid(True)
        raise NotImplementedError('TCF plotting not implemented yet')
        plt.show()
        sys.exit()
        return t,C


    def get_C_UDs(self):
        raise NotImplementedError('get_C_UDs not implemented yet')
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



