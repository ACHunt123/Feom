#!/usr/bin/env python3
import numpy as np
from Feom.src.initialize.setup import Setup,SimConfig
from pyA4.Bose_BCF import BoseBCF
from Feom.potentials.spin_boson import Spin_boson
from types import SimpleNamespace

# general params
L=4
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

### Setup the system
pot = Spin_boson(SimpleNamespace(ns=ns, Delta=Delta, eps=eps), None)

### Setup the bath
Xi_d = np.sqrt(Omega_D**2 - gamma_D**2)
Jw_pos_poles = [Xi_d + 1.j * gamma_D,  -Xi_d + 1.j * gamma_D   ]
res_val = (0.5* lambda_D * Omega_D**2) / (1.j * Xi_d)
Jw_pos_residues= [res_val, -res_val]




# generate the BCF
bcf = BoseBCF(beta=beta)
bcf.set_Jw(Jw_pos_poles, Jw_pos_residues)
C_ks,gam_ks,zeta = bcf.compute_mats_infin_bcf(K)
zeta=0


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
Vcross = np.kron(pot.s_mat,I) - np.kron(I,pot.s_mat.T)  # commutator superoperator for the system-bath coupling operator
Xi= -1 * (zeta/2) * Vcross @ Vcross # Add on the terminator contribution from the delta function in BCF

### Build the dictionaries
sys_args = {
    's_mat': pot.s_mat,
    'H_mat': pot.H_mat,}

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



