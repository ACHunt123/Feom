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
import copy

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
assert params.bathmode == 'AAA'
cothpoles = Setup(params)
cothpoles.params.mu=params.K # needed as params.mu is used in the cothAAA.py file

## COMMAND HERE TO GET THE AAA DECOMPOSITION OF THE COTH FUNCTION

def Sbeta(w,mode='exact',M=1):
   ''' Calculate Sbeta (fourier tranform of the BCF)
   either exactly, or using the AAA poles
   we have to treat the divergence of coth at w=0 carefully'''
   if mode=='exact':
      coth_term_times_w = 2/(params.beta*params.hbar) + (4*w**2/(params.beta*params.hbar))*(cothpoles.bath.P(w))
   elif mode=='AAA coth poles':
      coth_term_times_w = 2/(params.beta*params.hbar) + (4*w**2/(params.beta*params.hbar))*(cothpoles.bath.P_aaa_realcoeffs(w)+cothpoles.bath.k)
   elif mode=='M low freq approx':
      coth_term_times_w = 2/(params.beta*params.hbar) + (4*w**2/(params.beta*params.hbar))*(cothpoles.bath.P(w,M=M))
   else: sys.exit('invalid mode for Sbeta function')
   Jw_over_w = params.eta*params.gam/(w**2+params.gam**2)
   return (coth_term_times_w + w)*Jw_over_w

### Run the AAA decomposition directly on Sbeta (Xu et al)
Xumethod = copy.deepcopy(cothpoles)
Xumethod.bath.mu = cothpoles.bath.mu//2 + 1 # only half the number of poles can be used for Sbeta (they are complex)
_,support = Xumethod.bath.get_support_and_values(mode='log',N_support = 50000) # Get support used for the cothpoles decomposition
values= Sbeta(support)         
print('old mu = ',cothpoles.bath.mu) # store the old mu value (number of coth poles)
print('new mu = ',Xumethod.bath.mu)
_,_,Xumethod.bath.konstant,_=get_coeffs(Xumethod.bath,support,values,ext_fname='_Sbeta_decomp') 

## COMMAND HERE TO GET THE AAA DECOMPOSITION OF SBETA

def Sbeta_AAA(w):
   '''Calculate Sbeta using the poles from the AAA decomposition of the ENITRE FUNCTION'''
   result = 0.0+0j if np.isscalar(w) else np.zeros_like(w,dtype=np.complex128)
   w=np.array(w,dtype=np.complex128)
   for kk in range(len(Xumethod.bath.res_original)):
      result += Xumethod.bath.res_original[kk] / (w - Xumethod.bath.poles_original[kk])
   return result + Xumethod.bath.konstant


### Approximate Sbeta using the low frequency approx (finding the maximu, M needed to get a good approx)
# here we use the same setup as the initial coth decomposition``
cothpoles_M = copy.deepcopy(cothpoles)
max_dev= 1e-3*cothpoles_M.bath.gam # max deviation
M=1
while True:
   delta=np.max(np.abs(Sbeta(support,'M low freq approx',M) - Sbeta(support,'exact')))
   M+=5
   if delta>max_dev:
      M-=5
      break
values = cothpoles_M.bath.P(support,M=M)
_,_,cothpoles_M.bath.konstant,_=get_coeffs(cothpoles_M.bath,support,values,ext_fname=f'_M_{M}_low_freq_approx') 

## COMMAND HERE TO GET THE AAA DECOMPOSITION OF SBETA LOW FREQ APPROX
#  1i*imag(res(j)) ./ (xx - 1i*imag(pol(j)));
def Sbeta_M_aprox(w):
   '''Calculate Sbeta using the poles from the AAA decomposition of the pole function with the low frequency approx'''
   def P_AAA_M_approx(w):
        result = 0.0 if np.isscalar(w) else np.zeros_like(w,dtype=np.complex128)
        for k in range(len(cothpoles_M.bath.res_original)):
            result += 1j*np.imag(cothpoles_M.bath.res_original[k]) / (w - 1j*np.imag(cothpoles_M.bath.poles_original[k]))
        return result
        
   coth_term_times_w = 2/(params.beta*params.hbar) + (4*w**2/(params.beta*params.hbar))*(P_AAA_M_approx(w)+cothpoles_M.bath.konstant)
   Jw_over_w = params.eta*params.gam/(w**2+params.gam**2)
   return (coth_term_times_w + w)*Jw_over_w




fig, (ax) = plt.subplots(1,1,figsize=(8,10))

ax.plot(support,Sbeta(support,'AAA coth poles'),label=r'AAA coth (imaginary) poles K='+str(cothpoles.params.mu))
ax.plot(support,Sbeta_AAA(support).real,label=r'AAA complex poles N='+str(Xumethod.bath.mu))
ax.plot(support,Sbeta(support,'exact'),ls='--',color='k',label=r'Exact $S_{\beta}(\omega)$')
ax.plot(support,Sbeta_M_aprox(support).real,label=r'Low freq approx M='+str(M)+f' (difference is {delta:.2e})')

# ax1.plot(support,Sbeta(support,'M low freq approx',M),ls='-',label=r'$S_{\beta}(\omega) M=$'+str(M)+f' (difference is {delta:.2e})')
# ax1.plot(support,Sbeta(support,'exact'),ls='--',label=r'$S_{\beta}(\omega) $ exact')
# ax2.plot(support,(bath.P(support,M=M)),ls='-',label=r'$P(\omega)$ M='+str(M))
# ax2.plot(support,(bath.P(support)),ls='--',label=r'$P(\omega)$ exact')
# ax1.set_xlim(-100,100)
# ax2.set_xlim(-0.1,0.1)
# ax1.legend()
# ax2.legend()
# ax.plot(support,Sbeta(support,'exact'),ls='--',color='k',label=r'Exact $S_{\beta}(\omega)$')
ax.set_xlabel(r'Frequency $(\omega)$')
ax.set_ylabel(r'$S_{\beta}(\omega)$')
ax.set_title(r'$S_{\beta}(\omega)$'+f' for Debye Bath with beta={params.beta}.\n Setup such that both HEOMS have the same number of ADOs indices.')
ax.set_xlim(-100,100)

plt.legend()
#save teh figure
plt.tight_layout()
# plt.savefig(f'Sbeta_comparison_beta{params.beta}_K{oldmu}.pdf')
plt.show()
sys.exit()


# Pw[i] = (self.beta*self.hbar/4)*(1/wi)*(1/np.tanh(self.beta*self.hbar*wi/2) - 2/(self.beta*self.hbar*wi))



# run_aaa_fromfile(0,'/home/colehunt/data/test/aaa_K0','/home/colehunt/data/test/aaa_K0/aaa_data_quadrature_nw25001.txt','_quadrature_nw25001.txt',false)
# function run_aaa_fromfile(K, location,filename,extension,terminate)


plt.show()

# get_coeffs(params, support, bath.J(support))






