#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from Feom.src.fermions.initialize.setup import Setup, SimConfig
from Feom.src.fermions.hierarchy.fermion_ops import JW_rl_ops

# general params
L = np.inf
# K = 4
beta = 4
tmax = 100
dt = 0.001

# system params
n_sites = 1
ed_s = [1] # the energies for the sites

### Setup the system
c_dags, c_s, ns = JW_rl_ops(n_sites)
H_mat = np.zeros((ns, ns), dtype=complex)
for site_i in range(n_sites):
    H_mat += ed_s[site_i] * c_dags[site_i] @ c_s[site_i]

# perturbation matrices 
V_ks_plus = [c_dags[0] for k in range(2)]

### Setup the bath (single debye mode + 1 matsubara)
W = 1.0           # Bath bandwidth (cutoff frequency)
Gamma = 0.1       # Coupling strength
matsubara_1 = np.pi / beta  # approx 1.5708

# frequencies
gam_ks_plus = np.array([W, matsubara_1], dtype=complex)
gam_ks_mnus = np.array([W, matsubara_1], dtype=complex)

# prefactors (FIXED: f(w) vs 1-f(w) consistency)
fermi_at_minus_iW = 1.0 / (np.exp(-1j * beta * W) + 1.0)

c_W_plus = (Gamma * W / 2.0) * fermi_at_minus_iW
c_W_mnus = (Gamma * W / 2.0) * (1.0 - fermi_at_minus_iW)

c_M_plus = -1j * (Gamma * W**2) / (beta * (W**2 - matsubara_1**2))
c_M_mnus = -c_M_plus # Matsubara residue is inverted for the emission correlation

C_ks_plus = np.array([c_W_plus, c_M_plus], dtype=complex)
C_ks_mnus = np.array([c_W_mnus, c_M_mnus], dtype=complex)

### Build the dictionaries
sys_args = {
    'V_ks_plus': V_ks_plus,
    'H_mat': H_mat,
}

bath_args = {
    'C_ks_plus': C_ks_plus,
    'C_ks_mnus': C_ks_mnus,
    'gam_ks_plus': gam_ks_plus,
    'gam_ks_mnus': gam_ks_mnus,
}

params_args = {
    'dt': dt,
    'tmax': tmax,
    'L': L,
}

# Initialize
my_config = SimConfig(sys_args, bath_args, params_args)
sim = Setup(my_config)

# set initial conditions
rhos_0 = np.zeros((2, 2), dtype=complex)
rhos_0[1, 1] = 1.0 # Population in state 1
sim.set_initial_ADOs(rhos_0, '0th')

# Run
sim.generate_input_files()
sim.insert_executable()
sim.go(cleanup=1)

rho11 = sim.rho[:, 1, 1]  
rho00 = sim.rho[:, 0, 0]
rho10 = sim.rho[:, 1, 0]
rho01 = sim.rho[:, 0, 1]


# =====================================================================
# EXACT ANALYTICAL CALCULATION FOR VERIFICATION
# =====================================================================


def compute_analytical_rho11(t_arr, W, Gamma, ed, C_mnus, gam_mnus):
    # Roots of the self-energy propagator denominator: s^2 + b*s + c = 0
    b = W + 1j * ed
    c = 1j * ed * W + (Gamma * W) / 2.0
    s1 = (-b + np.sqrt(b**2 - 4*c)) / 2.0
    s2 = (-b - np.sqrt(b**2 - 4*c)) / 2.0
    
    # Residues for u(t)
    A1 = (s1 + W) / (s1 - s2)
    A2 = (s2 + W) / (s2 - s1)
    
    s = [s1, s2]
    A = [A1, A2]
    
    rho11_exact = []
    for t in t_arr:
        if t == 0:
            rho11_exact.append(1.0)
            continue
            
        # Homogeneous part (decay of initial occupation)
        u_t = A[0] * np.exp(s[0] * t) + A[1] * np.exp(s[1] * t)
        val = np.abs(u_t)**2
        
        # Inhomogeneous part (noise integration)
        noise = 0.0 + 0.j
        for j in range(2):
            for l in range(2):
                term_greater = 0.0 + 0.j
                term_less = 0.0 + 0.j
                
                sj_c = s[j].conj()
                sl = s[l]
                den_all = sj_c + sl
                exp_all = np.exp(den_all * t)
                
                for ck, gk in zip(C_mnus, gam_mnus):
                    # --- Region 1: tau2 > tau1 ---
                    den_sl_gk = sl - gk
                    den_sj_gk = sj_c + gk
                    exp_sl_gk = np.exp(den_sl_gk * t)
                    
                    # Defend against accidental division-by-zero limits
                    part1_g = (exp_all - exp_sl_gk) / den_sj_gk if np.abs(den_sj_gk) > 1e-9 else t * exp_all
                    part2_g = (exp_all - 1.0) / den_all if np.abs(den_all) > 1e-9 else t
                    
                    if np.abs(den_sl_gk) > 1e-9:
                        term_greater += (ck / den_sl_gk) * (part1_g - part2_g)
                    else:
                        term_greater += ck * (t**2 / 2.0)
                    
                    # --- Region 2: tau1 > tau2 ---
                    den_sj_gk_minus = sj_c - gk
                    den_sl_gk_plus = sl + gk
                    exp_sj_gk_minus = np.exp(den_sj_gk_minus * t)
                    
                    part1_l = (exp_all - exp_sj_gk_minus) / den_sl_gk_plus if np.abs(den_sl_gk_plus) > 1e-9 else t * exp_all
                    
                    if np.abs(den_sj_gk_minus) > 1e-9:
                        term_less += (np.conj(ck) / den_sj_gk_minus) * (part1_l - part2_g)
                    else:
                        term_less += np.conj(ck) * (t**2 / 2.0)
                        
                noise += np.conj(A[j]) * A[l] * (term_greater + term_less)
                
        val += noise.real
        rho11_exact.append(val)
        
    return np.array(rho11_exact)
rho11_analytic = compute_analytical_rho11(sim.t_arr, W, Gamma, ed_s[0], C_ks_mnus, gam_ks_mnus)

# =====================================================================
# PLOTTING AND COMPARISON
# =====================================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

axes[0].plot(sim.t_arr, rho11.real, 'o', color='crimson', label='Fermionic HEOM', markersize=4)
axes[0].plot(sim.t_arr, rho11_analytic, '-', color='navy', linewidth=2, label='Analytical Exact')
axes[0].set_xlabel('Time / a.u.', fontsize=12)
axes[0].set_ylabel('Site Population', fontsize=12)
axes[0].set_title('Site Population Validation', fontsize=13, fontweight='bold')
axes[0].grid(True, linestyle='--', alpha=0.6)
axes[0].legend(fontsize=11)

absolute_error = np.abs(rho11.real - rho11_analytic)
axes[1].semilogy(sim.t_arr, absolute_error, color='purple', linewidth=2, label='Abs Error')
axes[1].set_xlabel('Time / a.u.', fontsize=12)
axes[1].set_ylabel('Absolute Error', fontsize=12)
axes[1].set_title('HEOM Convergence Accuracy', fontsize=13, fontweight='bold')
axes[1].grid(True, linestyle='--', alpha=0.6)
axes[1].legend(fontsize=11)

plt.tight_layout()
plt.show()