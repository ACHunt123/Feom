import numpy as np
import matplotlib.pyplot as plt
# Default parameters
w_max_default=500
N_support_default=100000

from pyA4.Bose_BCF import BoseBCF
from pyA4.PyA4 import  A4Decomposition
def A4_BCF(Jw_pos_poles, Jw_pos_residues,W,beta,doplot=False,w_max=w_max_default,N_support=N_support_default):
    K = W-len(Jw_pos_poles) #A4 adds K poles to the BCF
    bcf = BoseBCF(beta=beta)
    # do the A4
    A4decomp=A4Decomposition(beta=beta,
                             K=K,
                             distribution='Bose',
                             N_support=N_support,
                             rational_decomposition_type='AAA',
                             w_max=w_max)
    eta_n,k_n = A4decomp.compute(doplot=doplot)
    # set the Rg and Jw poles/residues
    bcf.set_Jw(Jw_pos_poles, Jw_pos_residues)
    bcf.set_Rg_lorentzian_form(eta_n,k_n)
    # compute 
    C_ks,gam_ks,zeta = bcf.compute_bcf()
    return C_ks,gam_ks,zeta

from scipy.interpolate import AAA as scipy_AAA
def A3_BCF(Jw_pos_poles, Jw_pos_residues,W,beta,doplot=False,w_max=w_max_default,N_support=N_support_default):
    bcf = BoseBCF(beta=beta)
    bcf.set_Jw(Jw_pos_poles, Jw_pos_residues)
    Sbeta,omega= Sbeta_exact(Jw_pos_poles, 
                             Jw_pos_residues,
                             beta,
                             w_max=w_max,
                             N_support=N_support)
    r=scipy_AAA(omega, Sbeta, max_terms=W+1,rtol=0)
    poles,residues=r.poles(),r.residues()
    if doplot:
        plt.plot(omega,Sbeta,label='Sbeta exact')
        plt.plot(omega,r(omega),label='Sbeta AAA')
        plt.show()
    # fourier transform back to get the AAA BCF approximation
    mask=np.imag(poles)>0
    poles=poles[mask]
    residues=residues[mask]

    C_ks = 2*np.pi*1.j*residues[:]
    gam_ks = -1.j*poles[:]
    zeta=r(np.inf)
    return C_ks,gam_ks,zeta

def Sbeta_exact(Jw_pos_poles, Jw_pos_residues,beta,w_max=w_max_default,N_support=N_support_default):
    bcf = BoseBCF(beta=beta)
    bcf.set_Jw(Jw_pos_poles, Jw_pos_residues)
    omega=np.linspace(-w_max,w_max,N_support)+1e-15
    Sbeta=(1/np.pi)*(bcf.J(omega)*(1/np.tanh(beta*omega/2)/2 - (1/2)))
    return Sbeta,omega


def make_bcf_func(C_ks, gam_ks, zeta=0):
    # make the bcf and ignore the delta function
    def bcf(t):
        t_arr = np.atleast_1d(t)
        vals = np.sum(C_ks[:,None] * np.exp(-gam_ks[:,None] * t_arr[None,:]), axis=0)
        return vals if not np.isscalar(t) else vals[0]
    return bcf
    