#!/usr/bin/env python3
'''
   +---------------------------------------+
   |   FEOM: Fortran heirarchical          |  
   |       Equations Of Motion             |
   |           By A. C. Hunt 2025          |
   +---------------------------------------+
'''
from Feom.setup import Setup
from Feom.parser import params
import matplotlib.pyplot as plt
import sys,os
import numpy as np
from Feom.baths.coth_decomp.cothAAA import get_coeffs

###Store the command used to run this script in a file for later use
# Get the name of the Python interpreter (e.g., 'python' or 'python3')
python_cmd = os.path.basename(sys.executable)
# Reconstruct the command used to run the script
full_command = f"{python_cmd} {' '.join(sys.argv)}"
# Write it to a file named run.sh
with open("run.sh", "w") as f:
    f.write("#!/bin/bash\n")
    f.write(full_command + "\n")
####

### Load all the parameters into the setup object  (will perform the AAA decomposition if needed)
sim = Setup(params)
bath=sim.bath
params=sim.params
params.mu=params.K # needed as params.mu is used in the cothAAA.py file
# after this point, there will be a command to run AAA on the pole function

support=np.linspace(-100,100,1000)
fig,ax=plt.subplots()
def Sbeta(w,mode='exact',M=1):
   ''' Calculate Sbeta (fourier tranform of the BCF)
   either exactly, or using the AAA poles
   we have to treat the divergence of coth at w=0 carefully'''
   if mode=='exact':
      coth_term_times_w = 2/(params.beta*params.hbar) + (4*w**2/(params.beta*params.hbar))*(bath.P(w)+bath.k)
   elif mode=='AAA coth poles':
      coth_term_times_w = 2/(params.beta*params.hbar) + (4*w**2/(params.beta*params.hbar))*(bath.P_aaa_realcoeffs(w)+bath.k)
   elif mode=='M low freq approx':
      coth_term_times_w = 2/(params.beta*params.hbar) + (4*w**2/(params.beta*params.hbar))*(bath.P(w,M=M)+bath.k)
   else: sys.exit('invalid mode for Sbeta function')
   Jw_over_w = params.eta*params.gam/(w**2+params.gam**2)
   return (coth_term_times_w + w)*Jw_over_w

### Generate a command to run the AAA decomposition on Sbetaw
# Generate the support and values for the AAA decomposition (copy pasta)
N_support = 25000
eps=1e-4
w_max=bath.gam # start with the cuttoff frequency
Jw_min_tol = 1e-5      # tolerance for the maximum frequency of the grid for the AAA decomposition
while bath.J(w_max) > Jw_min_tol: w_max += 10 # find the maximum frequency where J(w) is still non-zero
x_pos = np.logspace(np.log10(eps), np.log10(w_max), N_support // 2)
support = np.concatenate((-x_pos[::-1], [0.0], x_pos))
values= Sbeta(support) # values of Sbeta function at the support points
oldmu=bath.mu
print('old mu = ',bath.mu)
bath.mu=bath.mu//2 + 1
print('new mu = ',bath.mu)
# get the AAA decomposition of Sbeta OVERWRITES THE POLES IN THE BATH OBJECT
_,_,konstant,_=get_coeffs(bath,support,values,ext_fname='_Sbeta_decomp') 
def Sbeta_AAA(w):
   '''Calculate Sbeta using the poles from the AAA decomposition of the ENITRE FUNCTION'''
   result = 0.0+0j if np.isscalar(w) else np.zeros_like(w,dtype=np.complex128)
   w=np.array(w,dtype=np.complex128)
   for kk in range(len(bath.res_original)):
      result += bath.res_original[kk] / (w - bath.poles_original[kk])

   return result + konstant



# ax.plot(support,Sbeta(support,'AAA coth poles'),label=r'AAA coth (imaginary) poles K='+str(oldmu))
# ax.plot(support,Sbeta_AAA(support).real,label=r'AAA complex poles N='+str(bath.mu))
# ax.plot(support,Sbeta(support,'exact'),ls='--',color='k',label=r'Exact $S_{\beta}(\omega)$')
fig, (ax1,ax2) = plt.subplots(2,1,figsize=(8,10))
for M in [10,50,100,200,500]:
   delta=np.max(np.abs(Sbeta(support,'M low freq approx',M) - Sbeta(support,'exact')))
   ax1.plot(support,Sbeta(support,'M low freq approx',M),ls='-',label=r'$S_{\beta}(\omega) M=$'+str(M)+f' (difference is {delta:.2e})')
   ax2.plot(support,(bath.P(support,M=M)),ls='-',label=r'$P(\omega)$ M='+str(M))
ax1.set_xlim(-100,100)
ax2.set_xlim(-0.1,0.1)
ax1.legend()
ax2.legend()
# ax.plot(support,Sbeta(support,'exact'),ls='--',color='k',label=r'Exact $S_{\beta}(\omega)$')
ax.set_xlabel(r'Frequency $(\omega)$')
ax.set_ylabel(r'$S_{\beta}(\omega)$')
ax.set_title(r'$S_{\beta}(\omega)$'+f' for Debye Bath with beta={params.beta}.\n Setup such that both HEOMS have the same number of ADOs indices.')
ax.set_xlim(-100,100)

plt.legend()
#save teh figure
plt.tight_layout()
plt.savefig(f'Sbeta_comparison_beta{params.beta}_K{oldmu}.pdf')
plt.show()
sys.exit()


# Pw[i] = (self.beta*self.hbar/4)*(1/wi)*(1/np.tanh(self.beta*self.hbar*wi/2) - 2/(self.beta*self.hbar*wi))



# run_aaa_fromfile(0,'/home/colehunt/data/test/aaa_K0','/home/colehunt/data/test/aaa_K0/aaa_data_quadrature_nw25001.txt','_quadrature_nw25001.txt',false)
# function run_aaa_fromfile(K, location,filename,extension,terminate)


plt.show()

# get_coeffs(params, support, bath.J(support))






