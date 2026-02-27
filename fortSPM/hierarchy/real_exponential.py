from scipy.sparse import csr_matrix
import numpy as np

def generate_liouvillian(sim):


    H= sim.pot.H_mat
    S= sim.pot.s_mat
    I=np.eye(sim.params.ns)

    Vcross = np.kron(S,I) - np.kron(I,S.T)  # commutator superoperator for the system-bath coupling operator
    Xi= -2 * (1/2) * Vcross @ Vcross # Add on the terminator contribution from the delta function in BCF

    
    A_dense = -1.j*(np.kron(H,I) - np.kron(I,H.T)) + Xi

    # Convert the masked dense matrix into a Scipy CSR matrix
    A = csr_matrix(A_dense, dtype=np.complex128)

    sim.Liouvillian=A
