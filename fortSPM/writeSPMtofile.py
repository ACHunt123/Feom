import numpy as np
from scipy.sparse import csr_matrix

# Example sparse, complex matrix
A = csr_matrix([
    [100+0.j, 1+20.j, 0+0.j],
    [0+0.j, 200+0.j, 0+0.j],
    [300+0.j, 0+0.j, 0+0.j]
], dtype=np.complex128)

# Use 1-based indexing for Fortran
row_ptr = A.indptr + 1
col_ind = A.indices + 1
values  = A.data

with open("csr_matrix_fortran.txt", "w") as f:
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

print('all done')