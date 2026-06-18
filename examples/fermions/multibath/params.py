import numpy as np
from Feom.src.fermions.hierarchy.fermion_ops import JW_rl_ops

# --- General Params ---
L = np.inf # Truncation tier
n_sites = 2
tmax = 1000
dt = 0.1
ed_s = [1.0, 1.5]  # Different energies for site 0 and site 1
t_hop = 0.5        # Hopping parameter (if you want the sites to talk to each other)

# --- System Setup ---
c_dags, c_s, ns = JW_rl_ops(n_sites)
H_mat = np.zeros((ns, ns), dtype=complex)

# On-site energies
for i in range(n_sites):
    H_mat += ed_s[i] * (c_dags[i] @ c_s[i])
# Hopping 
H_mat += t_hop * (c_dags[0] @ c_s[1] + c_dags[1] @ c_s[0])

# --- Perturbation Matrices --- (1 modes on each bath)
V_ks_plus = [c_dags[0],c_dags[1]]
V_ks_mnus = [c_s[0],c_s[1]]

# --- Bath Parameters ---
W     = [0.10, 0.10]    # Narrower Lorentzian = very long memory time (tau = 1/W = 10)
Omega = [1.0, 1.0]     # Shifting the bath peak to compete with system hopping
Gamma = [0.2, 0.2]     # Moderate coupling
mu    = [0.5, 0.5]     # Chemical potential near system energies
beta  = 1.0            #

C_ks_plus, gam_ks_plus = [], []
C_ks_mnus, gam_ks_mnus = [], []

for i in range(n_sites):
    # Plus Branch (c_dag)
    z_plus = Omega[i] - 1.0j * W[i]
    fermi_plus = 1.0 / (np.exp(beta * (z_plus - mu[i])) + 1.0)
    C_ks_plus.append(Gamma[i] * W[i] * (1.0 - fermi_plus))
    gam_ks_plus.append(W[i] + 1.0j * Omega[i])
    
    # Minus Branch (c)
    z_mnus = Omega[i] + 1.0j * W[i]
    fermi_mnus = 1.0 / (np.exp(beta * (z_mnus - mu[i])) + 1.0)
    C_ks_mnus.append(Gamma[i] * W[i] * fermi_mnus)
    gam_ks_mnus.append(W[i] - 1.0j * Omega[i])

# Convert to numpy arrays as expected by your Liouvillian generator
C_ks_plus = np.array(C_ks_plus, dtype=complex)
gam_ks_plus = np.array(gam_ks_plus, dtype=complex)
C_ks_mnus = np.array(C_ks_mnus, dtype=complex)
gam_ks_mnus = np.array(gam_ks_mnus, dtype=complex)