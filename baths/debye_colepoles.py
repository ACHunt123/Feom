import numpy as np
import matplotlib.pyplot as plt
import sys, os 
# Get the directory of the current script file
script_dir = os.path.dirname(os.path.abspath(__file__))

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
        self.aaa_tol = 1e-10        # tolerance for the AAA decomposition
        self.Jw_min_tol = 1e-2      # tolerance for the maximum frequency of the grid for the AAA decomposition
        self.minres_tol = 1e-6      # tolerance for the minimum abs value of a residue in the AAA decomposition

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
        self.get_coefs()
        self.TCF(plotme=True,ax=plt) # Calculate the TCF for the bath and plot it
        params.K = len(self.w_i) # Update the number of exponentials in the params object to match the number of poles found
        params.lowTCorr=True
        ### Calculate the coefficients C_U, c_D_LEFT, c_D_RIGHT for the bath (that are used in the FEOM code)
        self.get_C_UDs()

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
        ### Calculate the proposed extent of the support such that J(w) has decayed to 0
        w_max=self.gam # start with the cuttoff frequency
        while self.J(w_max) > self.Jw_min_tol: w_max += 10 # find the maximum frequency where J(w) is still non-zero
        print(f'Maximum frequency for the AAA decomposition: {w_max} (tolerance {self.Jw_min_tol})')
        # support = np.linspace(-w_max,w_max,int(200*w_max),dtype=np.complex128) # support for the AAA decomposition
        support = np.linspace(-250,250,20000,dtype=np.complex128) # support for the AAA decomposition
        values = self.P(support)                     # values of the pole function at the support points
        ### Use the AAA decomposition to get the coefficients
        folder = f'aaa_K{self.mu}'                   # folder to save the aaa files
        if not os.path.exists(folder): os.makedirs(folder)
        data = np.column_stack((support.real, values.real, values.imag))  
        print(f'Saving support and values to {folder}/aaa_data.txt')
        np.savetxt(f'{folder}/aaa_data.txt', data, header='Support Re[Values] Im[Values]', comments='')         # save the support and values to be read by matlab
        ### Run the AAA decomposition in MATLAB
        if(0):                                      # run the MATLAB script using the system command
            print('Running AAA decomposition in MATLAB...')
            os.system(f"matlab -batch 'cd {script_dir}/aaa;run_aaa_fromfile({self.mu},{os.getcwd()}/{folder})' > /dev/null 2>&1")    # run the matlab script to get the AAA coefficients
            print('AAA decomposition complete, loading results...')
        else:                                       # print out the command to run the MATLAB script if it has not already been run
            if not os.path.exists(f'{folder}/pol_real.txt'):
                print('run the following command in MATLAB to get the AAA coefficients:')
                print(f"\nrun_aaa_fromfile({self.mu},'{os.getcwd()}/{folder}')")
                ### append the command to a file for later use
                with open(f'{script_dir}/aaa/commands_to_run.m', 'a') as f:
                    if self.mu>0: f.write(f"run_aaa_fromfile({self.mu},'{os.getcwd()}/{folder}')\n")
                print(f"\n")
                sys.exit()
            else:
                print('AAA decomposition already done, loading results from files...')

        ### Load the aaa results
        repoles = np.loadtxt(f'{folder}/pol_real.txt')
        impoles = np.loadtxt(f'{folder}/pol_imag.txt')
        self.poles = repoles + 1.j * impoles
        reres = np.loadtxt(f'{folder}/res_real.txt')
        imres = np.loadtxt(f'{folder}/res_imag.txt')
        self.k = np.loadtxt(f'{folder}/k.txt')  # load the constant shift k
        self.res = reres + 1.j * imres
        self.res_original = self.res.copy()         # save the original residues for later use
        self.poles_original = self.poles.copy()     # save the original poles for later use
        ### Clean up the poles and residues
        mask= np.abs(self.res)> self.minres_tol
        self.res = self.res[mask]      # remove any tiny residues
        self.poles = self.poles[mask]  # remove the corresponding poles
        self.res = np.imag(self.res)*1.j            # remove the real parts, as by symmetry they should be zero
        self.poles = np.imag(self.poles)*1.j
        ### Calulate the real coefficients w_i and gamma_i from conjugate pairs of poles and residues
        upper_poles= []; upper_res = []
        for k in range(len(self.poles)):
            if np.imag(self.poles[k]) > 0:
                upper_poles.append(self.poles[k])
                upper_res.append(self.res[k])
        upper_poles = np.array(upper_poles,dtype=np.complex128) ; upper_res = np.array(upper_res,dtype=np.complex128)
        self.w_i = np.imag(upper_poles)                             # these are the new prequencies
        self.gam_i = -2*np.imag(upper_poles)*np.imag(upper_res)     # these are the new gammas
        ### Recalculate the number of exponentials (as the AAA algorithm might have changed the number of poles)
        if self.N_exp != len(self.poles):
            print(f'Total of frequencies {len(self.w_i)+self.N_nonmats} does not match number of exponentials proposed ({self.N_exp}), changing now.')
            print(f'With this new set, K={len(self.w_i)}.')
        self.N_exp = len(self.w_i)+self.N_nonmats  # update the number of exponentials
        if(1): # Plot the approximated function and J(w)
            plt.figure(figsize=(10,5))
            plt.plot(support, values.real, label='Original Function', color='blue')
            plt.plot(support, self.P_aaa_realcoeffs(support).real+self.k, label='AAA Approximation, imaginary poles/residues', color='green')
            # plt.plot(support, self.P_aaa(support).real, label='AAA Approximation', color='red')
            plt.plot(support, self.J(support), label='J(w)', color='orange')
            plt.xlabel('Support')
            plt.ylabel('Function Value')
            plt.title('AAA Approximation of the Pole Function')
            plt.legend()
            plt.grid()
            plt.show()
        if(0): #plot the poles
            plt.figure(figsize=(12, 6))
            pole_ax = plt.subplot(121)
            pole_ax.set_title('Poles')
            pole_ax.set_xlabel('Real Part')
            pole_ax.set_ylabel('Imaginary Part')
            pole_ax.grid(True)
            pole_ax.plot(repoles, impoles, 'o', label=f'AAA, N={len(repoles)}', color='blue')
            pole_ax.legend()
            pole_ax.set_xlim(-1500, 1500)
            plt.show()


    # Calculate the C_ks and gam_ks
    def get_coefs(self):
        self.calc_poles()  # Calculate the poles and residues from the coth function, and set the w_i and gam_i attributes

        d = np.zeros(self.N_exp)
        ### Calculate the 0th term
        d0sum = np.sum(self.gam_i/(self.gam**2 - self.w_i**2))
        d[0] = self.eta/self.beta + (2*self.eta*self.gam**2/self.beta)*d0sum - 2*self.eta*self.gam**2*self.k/self.beta

        ### Calculate the rest of the terms (matsubara terms)
        d[1:] = -(2*self.eta*self.gam/self.beta) * self.w_i[:]*self.gam_i[:]/(self.gam**2*np.ones_like(self.w_i) - self.w_i[:]**2)
        dI = -self.hbar*self.eta*self.gam/2

        # Calculate the C_ks
        C_ks = np.zeros(self.N_exp,dtype=complex)
        C_ks[0] = (d[0] + dI*1.j)
        C_ks[1:] = d[1:]

        #Calculate the gam_ks 
        gam_ks = np.zeros(self.N_exp,dtype=complex)
        gam_ks[0] = self.gam
        gam_ks[1:] = self.w_i[:]

       
        self.C_ks = C_ks
        self.gam_ks = gam_ks

        self.lowTcoef=-2*self.eta*self.gam*self.k/(self.beta*self.hbar**2) 
        self.C0hot=None

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
            print(mode)
            plt.show()
        return t,C

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



