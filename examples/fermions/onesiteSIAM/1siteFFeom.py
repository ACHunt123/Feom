#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from Feom.src.fermions.initialize.setup import Setup, SimConfig

from params import *

### Build the dictionaries
sys_args = {
    'V_ks_plus': V_ks_plus,
    'H_mat': H_mat,}

bath_args = {
    'C_ks_plus': C_ks_plus,
    'C_ks_mnus': C_ks_mnus,
    'gam_ks_plus': gam_ks_plus,
    'gam_ks_mnus': gam_ks_mnus,}

params_args = {
    'dt': dt,
    'tmax': tmax,
    'L': L,
    'noSIA': 1,         # Whether to use SIA or not (DONT IF SIMULATING THE WHOLE HIERARCHY)
}

# Initialize
my_config = SimConfig(sys_args, bath_args, params_args)
sim = Setup(my_config)

# set initial conditions
from Feom.src.fermions.hierarchy.fermion_ops import generate_state
psi0=generate_state(1,[1])
rhos_0 = np.outer(psi0.conj(),psi0)
sim.set_initial_ADOs(rhos_0, '0th')

# Run
sim.generate_input_files()
sim.insert_executable()
sim.go(cleanup=1)

rho11=np.array([np.trace(c_dags[0] @ c_s[0]@sim.rho[it,:,:]) for it in range(len(sim.t_arr))])
trace=np.array([np.trace(sim.rho[it,:,:]) for it in range(len(sim.t_arr))])



qutipdata=np.loadtxt('data/rho11qutip.dat')
t_arr_qutip=qutipdata[:,0]
rho11_qutip=qutipdata[:,1]
trace_qutip=qutipdata[:,2]
# =====================================================================
# PLOTTING AND COMPARISON
# =====================================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

axes[0].plot(sim.t_arr, rho11.real, 'o', color='crimson', label='Fermionic HEOM', markersize=4)
axes[0].plot(sim.t_arr, trace, 'o', color='green', label='Fermionic HEOM trace', markersize=4)
axes[0].plot(t_arr_qutip, rho11_qutip, '-', color='navy', linewidth=2, label='Qutip')
axes[0].plot(t_arr_qutip, trace_qutip, '--', color='black', linewidth=2, label='Qutip tr')
axes[0].set_xlabel('Time / a.u.', fontsize=12)
axes[0].set_ylabel('Site Population', fontsize=12)
axes[0].set_title('Site Population Validation', fontsize=13, fontweight='bold')
axes[0].grid(True, linestyle='--', alpha=0.6)
axes[0].legend(fontsize=11)

try:
    absolute_error = np.abs(rho11.real - rho11_qutip.real)
    # Divide by the max value (or just 1.0) to prevent division by zero as it decays
    percent_error = 100 * absolute_error / np.max(np.abs(rho11.real))

    axes[1].plot(sim.t_arr, percent_error, color='purple', linewidth=2, label='Percentage Error')
    axes[1].set_xlabel('Time / a.u.', fontsize=12)
    axes[1].set_ylabel('Percent Error (%)', fontsize=12)
    axes[1].set_title('HEOM Convergence Accuracy', fontsize=13, fontweight='bold')
    axes[1].grid(True, linestyle='--', alpha=0.6)
    axes[1].legend()
except:
    print('align the t axes please')
plt.tight_layout()
plt.show()