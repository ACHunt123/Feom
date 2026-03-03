#!/usr/bin/env python3
import numpy as np
from Feom.src.initialize.setup import Setup,SimConfig
from pyA4.Bose_BCF import BoseBCF
from pyA4.PyA4 import  A4Decomposition
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
# do the A4
A4decomp=A4Decomposition(beta=beta,K=K,distribution='Bose',N_support=10000,rational_decomposition_type='AAA')
eta_n,k_n = A4decomp.compute(doplot=True)
# set the Rg and Jw poles/residues
bcf.set_Jw(Jw_pos_poles, Jw_pos_residues)
bcf.set_Rg_lorentzian_form(eta_n,k_n)
# compute 
C_ks,gam_ks,zeta = bcf.compute_bcf()

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



# 5. Run
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


import matplotlib.pyplot as plt
# calculate the exact answer by exponentiating the Liouvillian

from scipy.sparse.linalg import expm_multiply

# 1. Flatten your 2x2 initial density matrix into a 4x1 vector
# (This matches the vectorization of your Liouvillian superoperator)
rho_vec_0 = rhos_0.flatten()

# 2. Extract time array properties to feed into expm_multiply
t_start = sim.t_arr[0]
t_stop = sim.t_arr[-1]
n_steps = len(sim.t_arr)

# 3. Calculate the exact time evolution!
# expm_multiply is incredibly efficient and returns an array of shape (n_steps, 4)
exact_dynamics_vec = expm_multiply(
    sim.Liouvillian[0:ns**2,0:ns**2], 
    rho_vec_0[0:ns**2], 
    start=t_start, 
    stop=t_stop, 
    num=n_steps, 
    endpoint=True
)

# 4. Reshape the resulting vectors back into 2x2 density matrices for each time step
ns = sim.params.ns  # Your system size (2)
exact_dynamics = exact_dynamics_vec.reshape(n_steps, ns, ns)


# plt.plot(sim.t_arr,rho11)
# plt.plot(sim.t_arr,exact_dynamics[:, 1, 1].real,'--')
# plt.ylim(0,1)
# plt.xlim(0,tmax)
# plt.show()


# np.savetxt('manual_constantk.out',processed_data.real,header=data_labels)


