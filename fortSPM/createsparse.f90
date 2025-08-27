program create_mkl_sparse
    use mkl_spblas
    implicit none

    integer :: n_rows, n_cols, nnz
    integer, allocatable :: row_ptr(:), col_ind(:)
    real(8), allocatable :: values(:)
    type(mkl_sparse_matrix_handle) :: Ahandle
    integer :: status
    integer, parameter :: indexing = 1  ! 1 for one-based (Fortran style)
    integer :: iounit

    ! Open and read your CSR matrix from file (replace with your file reading code)
    open(unit=10, file='csr_matrix_with_header.txt', status='old')
    read(10, *) n_rows, n_cols, nnz
    allocate(row_ptr(n_rows+1))
    allocate(col_ind(nnz))
    allocate(values(nnz))
    do iounit=1, n_rows+1
        read(10, *) row_ptr(iounit)
    end do
    do iounit=1, nnz
        read(10, *) col_ind(iounit)
    end do
    do iounit=1, nnz
        read(10, *) values(iounit)
    end do
    close(10)

    ! Create MKL sparse CSR matrix handle
    status = mkl_sparse_d_create_csr(Ahandle, indexing, n_rows, n_cols, row_ptr, row_ptr(2:), col_ind, values)
    if (status /= 0) then
        print *, "Error creating MKL sparse matrix, status =", status
        stop
    end if

    print *, "MKL sparse matrix created successfully!"

    ! Now you can use Ahandle for operations like mkl_sparse_d_mv, etc.

    ! When done, destroy handle
    call mkl_sparse_destroy(Ahandle)

end program create_mkl_sparse
