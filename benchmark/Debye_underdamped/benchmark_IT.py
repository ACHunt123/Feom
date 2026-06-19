#!/usr/bin/env python3
import numpy as np
from Feom.src.bosons.initialize.setup import Setup,SimConfig
from pyA4.Bose_BCF import BoseBCF

# general params
L=2
K=3
ns=2
beta=10
tmax= 10
dt=0.01
# system params
Delta=2
eps=1


# bath params
# Tom's paramaterization
Omega_UBO = 1.0 
gamma_UBO = 0.2 
lambda_UBO = 0.5 
# convert to mine 
Omega_D = Omega_UBO
gamma_D = gamma_UBO/2
lambda_D = lambda_UBO
Xi_D = np.sqrt(Omega_D**2 - gamma_D**2)

### Setup the system
H_mat = np.zeros((ns,ns),dtype=complex)
H_mat[0,0] = -eps
H_mat[1,1] = eps
H_mat[1,0] = Delta
H_mat[0,1] = Delta
# perturbation matrix
q_mat = np.zeros((ns,ns),dtype=complex)
q_mat[0,0] = -1
q_mat[1,1] = 1

if(0): # do it the old way
    ### Setup the bath
    Jw_pos_poles = [Xi_D + 1.j * gamma_D,  -Xi_D + 1.j * gamma_D   ]
    res_val = (0.5* lambda_D * Omega_D**2) / (1.j * Xi_D)
    Jw_pos_residues= [res_val, -res_val]

    # generate the BCF
    bcf = BoseBCF(beta=beta)
    bcf.set_Jw(Jw_pos_poles, Jw_pos_residues)
    C_ks,gam_ks,zeta = bcf.compute_mats_infin_bcf(K)
    print(C_ks)
    print(gam_ks)
    print(zeta)

elif(0):    ### new one with analytically contindued debye baths
    gam1_DL = 1.j*(Xi_D - 1.j*gamma_D)
    gam2_DL = 1.j*(-Xi_D - 1.j*gamma_D)

    eta1_DL =1.j*(lambda_D*Omega_D**2/(2*Xi_D*gam1_DL))
    eta2_DL =-1.j*(lambda_D*Omega_D**2/(2*Xi_D*gam2_DL))


    Jw_pos_residues = [eta1_DL*gam1_DL,eta2_DL*gam2_DL]
    Jw_pos_poles=[1.j*gam1_DL,1.j*gam2_DL]
    bcf = BoseBCF(beta=beta)
    bcf.set_Jw(Jw_pos_poles, Jw_pos_residues)
    C_ks,gam_ks,zeta = bcf.compute_mats_infin_bcf(K)


else:    ### TWO separate analytically contindued debye baths
    gam1_DL = 1.j*(Xi_D - 1.j*gamma_D)
    gam2_DL = 1.j*(-Xi_D - 1.j*gamma_D)

    eta1_DL =1.j*(lambda_D*Omega_D**2/(2*Xi_D*gam1_DL))
    eta2_DL =-1.j*(lambda_D*Omega_D**2/(2*Xi_D*gam2_DL))

    Jw_pos_residues1 = [eta1_DL*gam1_DL]
    Jw_pos_poles1=[1.j*gam1_DL]
    Jw_pos_residues2 = [eta2_DL*gam2_DL]
    Jw_pos_poles2=[1.j*gam2_DL]

    bcf1 = BoseBCF(beta=beta)
    bcf1.set_Jw(Jw_pos_poles1, Jw_pos_residues1)
    C1_ks,gam1_ks,zeta1 = bcf1.compute_mats_infin_bcf(K)

    bcf2 = BoseBCF(beta=beta)
    bcf2.set_Jw(Jw_pos_poles2, Jw_pos_residues2)
    C2_ks,gam2_ks,zeta2 = bcf2.compute_mats_infin_bcf(K)

    # calculate them from the others
    zeta=zeta1+zeta2
    C_ks=np.zeros((K+2),dtype=complex)
    gam_ks=np.zeros((K+2),dtype=complex)
    gam_ks[0]=gam1_ks[0]
    gam_ks[1]=gam2_ks[0]
    gam_ks[2:]=gam2_ks[1:]

    C_ks[0]=C1_ks[0]
    C_ks[1]=C2_ks[0]
    C_ks[2:]=C2_ks[1:]+C1_ks[1:]


if(0):#plot the J(w)
    import matplotlib.pyplot as plt
    w = np.linspace(-5,5,1000)
    plt.plot(w,bcf.J(w),label='from poles and residues')
    def J_analytical(w):
        numerator = 4 * lambda_D * gamma_D * (Omega_D**2) * w
        denominator = (w**2 - Omega_D**2)**2 + 4*(gamma_D**2 * w**2)
        return numerator / denominator
    plt.plot(w,J_analytical(w),label='from coeffs',linestyle='--')
    plt.legend()
    plt.show()


    print('C_ks',C_ks)
    print('gam_ks',gam_ks)
    print(zeta)
    exit()


### Setup the terminator
I = np.eye(ns)
Vcross = np.kron(q_mat,I) - np.kron(I,q_mat.T)  # commutator superoperator for the system-bath coupling operator
Xi= -1 * (zeta/2) * Vcross @ Vcross # Add on the terminator contribution from the delta function in BCF
# Xi/=16
# LOOKS LIKE MY TERMINATOR IS 16 times LARGER THASNK TOMS

### Build the dictionaries
sys_args = {
    's_mat': q_mat,
    'H_mat': H_mat,}

bath_args = {
    'C_ks':C_ks,
    'gam_ks':gam_ks,
    'zeta':zeta}

params_args = {
    'dt': dt,
    'tmax': tmax,
    'L': L,
}

terminator_args = {
    'correction_type': 'same_for_each_ADO',
    'Xi':Xi}

# 4. Initialize
my_config=SimConfig(sys_args, bath_args, params_args,terminator_args)
sim = Setup(my_config)
# my_config.save('config.json')
# my_config=SimConfig.load('config.json')

# set initial conditions
rhos_0 = np.zeros((2, 2), dtype=complex)
rhos_0[1, 1] = 1.0 # Population in state 1
sim.set_initial_ADOs(rhos_0,'0th')

# Run
sim.generate_input_files()
sim.insert_executable()
sim.go(cleanup=1)

rho11 = sim.rho[:,1,1]  
rho00 = sim.rho[:,0,0]
rho10 = sim.rho[:,1,0]
rho01 = sim.rho[:,0,1]

processed_data = np.zeros((len(sim.t_arr),5),dtype=complex)
processed_data[:,0] = sim.t_arr
processed_data[:,1] = rho11 - rho00  # <s_z>
processed_data[:,2] = 1.j*(rho10 - rho01)  # <s_y>
processed_data[:,3] = (rho10 + rho01)  # <s_x>   
processed_data[:,4] = rho11 # Site 1 population
data_labels = '\n,Time /a.u. <s_z> <s_y> <s_x> P11'

np.savetxt('IT_test.out',processed_data.real,header=data_labels)



