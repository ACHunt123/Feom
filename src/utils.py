import numpy as np
import os
from scipy.sparse import csr_matrix

# Some useful functions for the HEOM code
# including reade-write functions and naming for FORTRAN executables

def write_sparse(filename,A):
    # Use 1-based indexing for Fortran
    row_ptr = A.indptr + 1
    col_ind = A.indices + 1
    values  = A.data
    with open(filename, "w") as f:
        # Write header, giving the size of the matrix and #nonzeros
        f.write(f"{A.shape[0]} {A.shape[1]} {A.nnz}\n")
        # Write row pointer: n_rows + 1 values
        for val in row_ptr:
            f.write(f"{val:d}\n")
        # Write column indices: nnz values
        for val in col_ind:
            f.write(f"{val:d}\n")
        # Write values with fixed formatting
        for val in values:
            f.write(f"{val.real:22.15e}\n".replace('e','d'))
            f.write(f"{val.imag:22.15e}\n".replace('e','d'))

def write_Zvec(filename, matrix):
    """
    Writes a numpy array to a text file for Fortran to read.
    Alternates real and imaginary parts using D22.15 format.
    """
    # Ensure the input is a complex numpy array
    matrix = np.asarray(matrix, dtype=np.complex128)
    # Create the header string containing the original shape
    # If 1D, it writes e.g., "4". If 2x2, it writes "2 2".
    shape_str = " ".join(map(str, matrix.shape))
    # Flatten the array for the 1D Fortran loop.
    # order='F' prevents transposition when read back into Fortran memory arrays.
    vec = matrix.flatten(order='F')
    with open(filename, "w") as f:
        f.write(f"{shape_str}\n") # Write the header (which Fortran read(70, *) will skip)
        f.write(f"{len(vec)}\n") # Write the number of elements
        # Write the alternating real and imaginary parts
        for val in vec:
            # Format nicely to match Fortran's (D22.15) requirements
            # standard python 'e' is replaced with Fortran 'D'
            real_str = f"{val.real:22.15e}".replace('e', 'D')
            imag_str = f"{val.imag:22.15e}".replace('e', 'D')
            
            f.write(f"{real_str}\n")
            f.write(f"{imag_str}\n")

def writeParams(filename,sim): # Writes small parameters into file
    params=sim.params
    bath=sim.bath
    nttot =int(params.tmax/params.dt)+1 # calculate the total number of time steps
    if not os.path.exists(f"tmp/"): os.makedirs(f"tmp/")
    with open(f"{filename}", "w") as f:
        f.write("ns,dt,nttot\n")
        f.write(f"{params.ns:10d}{params.dt:22.15e}{nttot:10d}\n".replace('e','d'))
        f.write("/\n")
    return

def FORT_SWITCHES(sim):
    params = sim.params
    ''' Supported compile-time switches (SWITCHES):
       -DPrint_ADOs         Print the ADOs to file every N timesteps
       -DSIA                Use SIA step instead of RK4 step (default)
     '''
    switches = []
    if params.print_ADOs:
        switches.append('Print_ADOs')
    if not params.noSIA:
        switches.append('SIA')
    # sort the switches to be alphabetical
    switches.sort()
    makefile_command= f'SWITCHES=" -D{" -D".join(switches)}"'
    executable_suffix='_'+'_'.join(switches)
    if len(switches) == 0:
        executable_suffix = ''
        makefile_command = 'SWITCHES=""'    
    return makefile_command, executable_suffix


if __name__=='__main__':
    

        # Define system size (2x2 density matrix means a 4x4 superoperator)
    N_sys = 2
    N_liouv = N_sys**2

    # ==========================================
    # 1. Random Complex Density Matrix (rho)
    # ==========================================
    # Generate random real and imaginary parts between 0.0 and 1.0
    rho_real = np.random.rand(N_sys, N_sys)
    rho_imag = np.random.rand(N_sys, N_sys)
    rho_matrix = rho_real + 1j * rho_imag

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


    ## write them all to files
    folder='matrices'
    write_sparse(f"{folder}/csr_matrix_fortran.txt",A)
    write_Zvec(f"{folder}/Fortrho_matrix.txt", rho_matrix)
    print('matrices written, now will do the product')

    # Flatten rho exactly as the text file writer does (Column-major / Fortran order)
    rho_vec = rho_matrix.flatten(order='F')

    # Perform the matrix-vector multiplication
    result_vec = A.dot(rho_vec)

    print("==========================================")
    print(" Python Result Vector (A * rho):")
    for i, val in enumerate(result_vec):
        print(f" [{i+1}] {val.real:15.5f}  + {val.imag:15.5f}j")
    print("==========================================")

