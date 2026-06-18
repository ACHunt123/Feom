import numpy as np
from Feom.src.fermions.hierarchy.fermion_ops import JW_rl_ops

# general params
L = np.inf
L = 2
tmax = 100
dt = 0.01

# system params
n_sites = 1
ed_s = [1] # the energies for the sites

### Setup the system
c_dags, c_s, ns = JW_rl_ops(n_sites)
H_mat = np.zeros((ns, ns), dtype=complex)
for site_i in range(n_sites):
    H_mat += ed_s[site_i] * c_dags[site_i] @ c_s[site_i]

# perturbation matrices 
V_ks_plus = [c_dags[0] for k in range(1)]
V_ks_mnus = [c_s[0] for k in range(1)]

import numpy as np

# --- Physical Parameters ---
W = 0.15        # Width of the bath Lorentzian peak (decay rate)
Omega = 2.0    # Center of the bath peak 
Gamma = 1.0    # Overall coupling strength
mu = 2.0       # Chemical potential 
beta = 0.01     # Inverse temperature (1/kT)

z_plus = Omega - 1.0j * W
fermi_plus = 1.0 / (np.exp(beta * (z_plus - mu)) + 1.0)

# C_plus uses 1 - f(z)
C_ks_plus = np.array([Gamma * W * (1.0 - fermi_plus)], dtype=complex)
gam_ks_plus = np.array([W + 1.0j * Omega], dtype=complex) 

# --- 2. The Minus Branch (s=-, lowering operator c) ---
# Closes in the Upper Half-Plane. Depends on bath ELECTRONS (f(z)).
z_mnus = Omega + 1.0j * W
fermi_mnus = 1.0 / (np.exp(beta * (z_mnus - mu)) + 1.0)

# C_mnus uses f(z)
C_ks_mnus = np.array([Gamma * W * fermi_mnus], dtype=complex)
gam_ks_mnus = np.array([W - 1.0j * Omega], dtype=complex)



