import numpy as np
from Feom.src.bosons.initialize.setup import Setup,SimConfig
from pyA4.Bose_BCF import BoseBCF
from scipy.linalg import expm

### Argparsed params
import argparse
def parse_args():
    parser = argparse.ArgumentParser(description="convergence params")
    # Convergence
    parser.add_argument("--dt", type=float, default=0.01)
    parser.add_argument("--teqm", type=float, default=500)
    parser.add_argument("--L", type=int, default=4)
    parser.add_argument("--K", type=int, default=3)
    parser.add_argument("--beta", type=float, default=10)
    parser.add_argument("--gam_DL", type=float, default=10)
    parser.add_argument("--ns", type=int, default=10)
    return parser.parse_args()
args=parse_args()
dt =  args.dt
L = args.L      
K = args.K      
beta = args.beta
teqm = args.teqm
gam_DL=args.gam_DL
ns=args.ns


### hardcoded Parameters
# general
tmax= 20
# bath
lambda_DL=0.05
eta_DL=lambda_DL*gam_DL
# system 
omega=1.



### Setup the system in energy eigenbasis
E_ns=(np.arange(ns)+0.5)*omega 
s_mat = np.zeros((ns,ns),dtype=complex)   #position basis
p_mat = np.zeros((ns, ns), dtype=complex) #momentum basis
def delta(i,j): return 1 if i==j else 0                 
for m in range(0,ns):
    for n in range(0,ns):
        s_mat[m,n] = np.sqrt(1/(2*omega))*(np.sqrt(n+1)*delta(m,n+1)+np.sqrt(n)*delta(m,n-1))
        p_mat[m,n] = 1j * np.sqrt(omega / 2) * (np.sqrt(n+1) * delta(m, n+1) - np.sqrt(n) * delta(m, n-1))
H0 = np.diag(E_ns)
# p_mat = 1.0j * (H0 @ s_mat - s_mat @ H0)
Hren = eta_DL*s_mat@s_mat/2.
H_mat = H0 + Hren #renormalization added

### Setup the bath
from pyA4.Bose_BCF import BoseBCF
from pyA4.PyA4 import  A4Decomposition
Jw_pos_residues = [eta_DL*gam_DL/2]
Jw_pos_poles=[1.j*gam_DL]
bcf = BoseBCF(beta=beta)
# do the A4
A4decomp=A4Decomposition(beta=beta,K=K,distribution='Bose',rational_decomposition_type='AAA',)
eta_n,k_n = A4decomp.compute(doplot=0)
# set the Rg and Jw poles/residues
bcf.set_Jw(Jw_pos_poles, Jw_pos_residues)
bcf.set_Rg_lorentzian_form(eta_n,k_n)
C_ks,gam_ks,zeta = bcf.compute_bcf()

### Setup the terminator REMOVED!! as it will diverge for white baths 
I = np.eye(ns)
Vcross = np.kron(s_mat,I) - np.kron(I,s_mat.T)  # commutator superoperator for the system-bath coupling operator
Xi= -1 * (zeta/2) * Vcross @ Vcross # Add on the terminator contribution from the delta function in BCF


### Build the dictionaries 
sys_args = {'s_mat': s_mat,'H_mat': H_mat,}
bath_args = {'C_ks':C_ks,'gam_ks':gam_ks,'zeta':zeta}
params_args_eqm = {'dt': dt,'tmax': teqm,'L': L}
params_args = {'dt': dt,'tmax': tmax,'L': L}
terminator_args = {'correction_type': 'same_for_each_ADO','Xi':Xi}


### Run the equilibriation
eqm_config=SimConfig(sys_args, bath_args, params_args_eqm,terminator_args)
sim = Setup(eqm_config)
# set initial conditions
rhos_0 = expm(-beta*H0)+0.j
rhos_0 /= np.trace(rhos_0)
sim.set_initial_ADOs(rhos_0,'0th')
sim.generate_input_files(tmp_folder='equilibriation')
sim.insert_executable()
sim.go(cleanup=0)
# record the populations for eqm analysis
processed_data = np.zeros((len(sim.t_arr),ns+1),dtype=complex)
processed_data[:,0] = sim.t_arr
for it in range(len(sim.t_arr)):
    processed_data[it,1:] = np.diagonal(sim.rho[it,:,:])
data_labels = '\n,Time /a.u. Populations 0->ns-1'
np.savetxt('equilibriation_populations.out',processed_data.real,header=data_labels)
# Read in the equilibriated ADOs and clean up
eqmADOs,rhos_eqm=sim.read_final_ADOs(include_rho_s=True)
np.save('rho_sys_eqm.npy', rhos_eqm)
sim.safe_cleanup()

### RUn the TCF
# perturb with the initial conditions rho_eq Q
config=SimConfig(sys_args, bath_args, params_args,terminator_args)
sim = Setup(config)
eqmADOs = eqmADOs.reshape(sim.params.Nados, ns, ns)
ADOs = np.einsum('ij,Kjk->Kik', s_mat, eqmADOs).flatten()
sim.set_initial_ADOs(ADOs,'all')
sim.generate_input_files(tmp_folder='after_equilibriation')
sim.insert_executable()
sim.go(cleanup=1)


processed_data = np.zeros((len(sim.t_arr),5),dtype=complex)
processed_data = np.zeros((len(sim.t_arr),ns+2),dtype=complex)
processed_data[:,0] = sim.t_arr
processed_data[:, 1] = np.einsum('ij, tji -> t', s_mat, sim.rho)
for it in range(len(sim.t_arr)):
    processed_data[it,2:] = np.diagonal(sim.rho[it,:,:])
data_labels = '\n,Time /a.u. Cqq Populations 0->ns-1'

np.savetxt('qtcf.out',processed_data.real,header=data_labels)
