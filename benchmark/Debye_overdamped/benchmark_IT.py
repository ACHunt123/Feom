#!/usr/bin/env python3
import numpy as np
from Feom.src.initialize.setup import Setup,SimConfig
from pyA4.Bose_BCF import BoseBCF

# general params
L=4
K=3
ns=2
beta=15
tmax= 10
dt=0.01
# system params
Delta=2
eps=1
# bath params
eta_DL=2
gam_DL=1


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


### Setup the bath
Jw_pos_residues = [eta_DL*gam_DL/2]
Jw_pos_poles=[1.j*gam_DL]
bcf = BoseBCF(beta=beta)
bcf.set_Jw(Jw_pos_poles, Jw_pos_residues)
C_ks,gam_ks,zeta = bcf.compute_mats_infin_bcf(K)

### Setup the terminator
I = np.eye(ns)
Vcross = np.kron(q_mat,I) - np.kron(I,q_mat.T)  # commutator superoperator for the system-bath coupling operator
Xi= -1 * (zeta/2) * Vcross @ Vcross # Add on the terminator contribution from the delta function in BCF

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



