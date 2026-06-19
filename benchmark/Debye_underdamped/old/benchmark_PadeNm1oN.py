#!/usr/bin/env python3
import numpy as np
from Feom.src.bosons.initialize.setup import Setup,SimConfig
from pyA4.Bose_BCF import BoseBCF
from pyA4.Pade import padeNm1oN
from Feom.potentials.spin_boson import Spin_boson
from types import SimpleNamespace


#  python3 
#  --L 3 --K 3 --ns 2 --beta 15
#  --hbar 1 --tmax 10.0 --dt 0.01 
# --potname spinboson --Delta 2.0 --eps 1.0 
# --bathname debyeCothpoles --bathmode Pade[N-1/N]
#  --eta 2 --gam 1 > /dev/null 
# general params
L=3
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
pot = Spin_boson(SimpleNamespace(ns=ns, Delta=Delta, eps=eps), None)

### Setup the bath
Jw_pos_residues = [eta_DL*gam_DL/2]
Jw_pos_poles=[1.j*gam_DL]
bcf = BoseBCF(beta=beta)
eta_n,k_n = padeNm1oN(K,beta)
bcf.set_Rg_lorentzian_form(eta_n,k_n)
bcf.set_Jw(Jw_pos_poles, Jw_pos_residues)

C_ks,gam_ks,zeta = bcf.compute_bcf()

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

np.savetxt('PadeNm1oN_test.out',processed_data.real,header=data_labels)



