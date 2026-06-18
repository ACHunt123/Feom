# main.py
import matplotlib.pyplot as plt
import qutip as qt
from qutip.solver.heom import HEOMSolver, FermionicBath

# Import all your variables and numpy matrices from params.py
from params import *

# QuTiP requires its own Qobj wrapper around numpy arrays to track Hilbert space dimensions
H_sys = qt.Qobj(H_mat)

# Define annihilation operators for BOTH sites
c_op_0 = qt.Qobj(c_s[0])
c_op_1 = qt.Qobj(c_s[1])
c_dag_op_0 = qt.Qobj(c_dags[0])

# construct the baths
bath0 = FermionicBath(
    c_op_0, 
    [C_ks_plus[0]],   [gam_ks_plus[0]], 
    [C_ks_mnus[0]],   [gam_ks_mnus[0]])

bath1 = FermionicBath(
    c_op_1, 
    [C_ks_plus[1]],   [gam_ks_plus[1]], 
    [C_ks_mnus[1]],   [gam_ks_mnus[1]])
# combine bath objs
bath = [bath0, bath1]
# Time array
t_arr = np.arange(0, tmax+dt, dt)

# start in |1><1|
rho0 = qt.basis(ns, 1) * qt.basis(ns, 1).dag()

# run it
solver = HEOMSolver(H_sys, bath, max_depth=L)

print(f"Running QuTiP HEOM for {len(t_arr)} time steps...")
# We track the population expectation value: <c^\dagger c>
I_sys = qt.qeye(ns) 

# Run the solver , tracking |0><0| and trace
result = solver.run(rho0, t_arr, e_ops=[c_dag_op_0 * c_op_0, I_sys])

rho11_qutip = result.expect[0]           # The population
qutip_traces = np.real(result.expect[1]) # The trace


print("Simulation complete.")

# save results for comparision
rho11_qutip = result.expect[0]
np.savetxt(
    'rho11qutip.dat', 
    np.column_stack((t_arr, rho11_qutip.real,qutip_traces)), 
    fmt='%.8e',               # Formats numbers in scientific notation
    header='Time Population', # Adds a nice header to the top of the file
    comments='# '             # Ensures the header is commented out
)
