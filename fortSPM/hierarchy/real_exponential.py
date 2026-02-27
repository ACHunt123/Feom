from scipy.sparse import csr_matrix
import numpy as np

def generate_liouvillian(sim):
    
    N_sys=2
    N_liouv = N_sys**2

    # for now just put in a placeholder
    # ==========================================
    # 2. Random Complex Sparse Matrix (A)
    # ==========================================
    # Generate a 4x4 dense random complex matrix first
    A_real = np.random.rand(N_liouv, N_liouv)
    A_imag = np.random.rand(N_liouv, N_liouv)
    A_dense = A_real + 1j * A_imag

    # Introduce sparsity by randomly setting ~50% of the elements to exactly zero
    # (This step makes it a "true" sparse matrix for testing)
    sparsity_mask = np.random.rand(N_liouv, N_liouv) > 0.5
    A_dense[sparsity_mask] = 0.0 + 0.0j

    # Convert the masked dense matrix into a Scipy CSR matrix
    A = csr_matrix(A_dense, dtype=np.complex128)

    sim.Liouvillian=A
