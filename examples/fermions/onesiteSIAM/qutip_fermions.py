# main.py
import matplotlib.pyplot as plt
import qutip as qt
from qutip.solver.heom import HEOMSolver, FermionicBath

# Import all your variables and numpy matrices from params.py
from params import *

# ==========================================
# 1. Convert Numpy Arrays to QuTiP Qobjs
# ==========================================
# QuTiP requires its own Qobj wrapper around numpy arrays to track Hilbert space dimensions
H_sys = qt.Qobj(H_mat)
c_op = qt.Qobj(c_s[0])
c_dag_op = qt.Qobj(c_dags[0])

# ==========================================
# 2. Construct Bath and Initial State
# ==========================================
# Pass the wrapped annihilation operator to the FermionicBath
bath = FermionicBath(c_op, C_ks_plus, gam_ks_plus, C_ks_mnus, gam_ks_mnus)

# Time array
t_arr = np.arange(0, tmax+dt, dt)

# Initial state: Assuming ns=2 (from your JW_rl_ops), we start with the site occupied.
# qt.basis(2, 1) creates the |1> state.
rho0 = qt.basis(ns, 1) * qt.basis(ns, 1).dag()

# ==========================================
# 3. Run HEOM Solver
# ==========================================
solver = HEOMSolver(H_sys, bath, max_depth=10)

print(f"Running QuTiP HEOM for {len(t_arr)} time steps...")
# We track the population expectation value: <c^\dagger c>
I_sys = qt.qeye(2) 

# 2. Run the solver with both operators in e_ops
result = solver.run(rho0, t_arr, e_ops=[c_dag_op * c_op, I_sys])


rho11_qutip = result.expect[0]           # The population
qutip_traces = np.real(result.expect[1]) # The trace



print("Simulation complete.")

# ==========================================
# 4. Extract and Plot Results
# ==========================================
rho11_qutip = result.expect[0]
np.savetxt(
    'data/rho11qutip.dat', 
    np.column_stack((t_arr, rho11_qutip.real,qutip_traces)), 
    fmt='%.8e',               # Formats numbers in scientific notation
    header='Time Population', # Adds a nice header to the top of the file
    comments='# '             # Ensures the header is commented out
)

# plt.figure(figsize=(8, 5))
# plt.plot(t_arr, rho11_qutip.real, label=r'$\rho_{11}$ (QuTiP HEOM)', color='blue', linewidth=2)
# plt.title('Fermionic HEOM Dynamics', fontsize=14)
# plt.xlabel('Time', fontsize=12)
# plt.ylabel('Population', fontsize=12)

# # Highly recommend keeping xlim tight to see the actual decay curve, 
# # otherwise 1000 time units will just look like a flat line at the steady state!
# plt.xlim(0, 50) 

# plt.grid(True, linestyle='--', alpha=0.7)
# plt.legend()
# plt.tight_layout()
# plt.show()