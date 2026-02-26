program test_main
    use complex_sparse_linalg
    use input_output
    implicit none

    type(complex_csr_matrix) :: A
    complex(8), allocatable :: vec(:),y(:)
    integer :: i
    ! Read in sparse matrix and dense vector (both complex(8))
    call read_matrix('matrices/csr_matrix_fortran.txt', A)
    call read_Zvec("matrices/Fortrho_matrix.txt", vec)
    ! Check all the dimensions etc.
    print *, "Checking dimensions..."
    if (size(vec) /= A%n_cols) then
        print *, "=========================================="
        print *, " FATAL ERROR: Dimension mismatch!"
        print *, " Matrix columns: ", A%n_cols
        print *, " Vector size:    ", size(vec)
        print *, "=========================================="
        stop "Execution aborted."
    else
        print *, " SUCCESS: Matrix and vector dimensions match! Size = ", size(vec)
    end if

    ! ! 2. Print the matrix to verify the contents
    ! call print_matrix(A)

    ! 3. Allocate output array
    allocate(y(size(vec)))

    ! 4. Perform the Matrix-Vector Multiplication!
    call csr_matvec_complex(A, vec, y)

    ! 5. Print the results to compare with Python
    print *, "=========================================="
    print *, " Fortran Result Vector (A * rho):"
    do i = 1, size(y)
        print '(" [",I0,"] ", F15.5, "  + ", F15.5, "j")', i, real(y(i)), aimag(y(i))
    end do
    print *, "=========================================="

    ! 3. Free the memory
    call destroy_matrix(A)
    deallocate(vec)

end program test_main


