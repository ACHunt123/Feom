module complex_sparse_linalg
    implicit none

    ! Define our native complex CSR matrix type
    type :: complex_csr_matrix
        integer :: n_rows
        integer :: n_cols
        integer :: nnz
        integer, allocatable :: row_ptr(:)
        integer, allocatable :: col_ind(:)
        complex(8), allocatable :: values(:)
    end type complex_csr_matrix

contains

    ! Subroutine to read the Python-generated file
    subroutine read_matrix(filename, A)
        character(len=*), intent(in)          :: filename
        type(complex_csr_matrix), intent(out) :: A
        integer :: i
        real(8) :: z_real, z_imag
        ! Open the file
        open(unit=10, file=filename, status='old', action='read')
        ! Read the header (n_rows, n_cols, nnz)
        read(10, *) A%n_rows, A%n_cols, A%nnz
        ! Allocate the CSR arrays based on the header
        allocate(A%row_ptr(A%n_rows + 1))
        allocate(A%col_ind(A%nnz))
        allocate(A%values(A%nnz))
        ! Read the row pointers
        do i = 1, A%n_rows + 1
            read(10, *) A%row_ptr(i)
        end do
        ! Read the column indices
        do i = 1, A%nnz
            read(10, *) A%col_ind(i)
        end do
        ! Read the complex values (real and imag on alternating lines)
        do i = 1, A%nnz
            read(10, '(D22.15)') z_real
            read(10, '(D22.15)') z_imag
            A%values(i) = dcmplx(z_real, z_imag)
        end do
        ! Close the file
        close(10)
        
    end subroutine read_matrix


   subroutine csr_matvec_complex(A, x, y)
        type(complex_csr_matrix), intent(in) :: A
        complex(8), intent(in)               :: x(:)  ! Changed to complex(8)
        complex(8), intent(out)              :: y(:)  ! Changed to complex(8)
        integer                              :: i, j
        complex(8)                           :: row_sum ! Changed to complex(8)
        
        ! Initialization is fast, but we can still parallelize it
        !$omp parallel workshare
        y = dcmplx(0.0d0, 0.0d0)  ! Initialize to complex zero
        !$omp end parallel workshare
        
        ! Parallelize over the rows. 
        ! 'guided' scheduling prevents load imbalance if row lengths vary wildly.
        !$omp parallel do private(i, j, row_sum) schedule(guided)
        do i = 1, A%n_rows
            
            ! Use a scalar to keep the running total in a fast CPU register
            row_sum = dcmplx(0.0d0, 0.0d0)
            
            ! Tell the compiler it is safe to vectorize this loop
            !DIR$ IVDEP 
            do j = A%row_ptr(i), A%row_ptr(i+1) - 1
                row_sum = row_sum + A%values(j) * x(A%col_ind(j))
            end do
            
            ! Write back to memory exactly once per row
            y(i) = row_sum
            
        end do
        !$omp end parallel do
        
    ! pretty print the matrix
    end subroutine csr_matvec_complex
    subroutine print_matrix(A)
        type(complex_csr_matrix), intent(in) :: A
        integer :: i, j
        
        print *, "=========================================="
        print *, "Sparse Matrix Details:"
        print *, "Dimensions: ", A%n_rows, "x", A%n_cols
        print *, "Non-zeros:  ", A%nnz
        print *, "------------------------------------------"
        do i = 1, A%n_rows
            do j = A%row_ptr(i), A%row_ptr(i+1) - 1
                print '(" Row: ", I0, "  Col: ", I0, "  Val: (", ES12.4, " , ", ES12.4, " i)")', &
                      i, A%col_ind(j), real(A%values(j)), aimag(A%values(j))
            end do
        end do
        print *, "=========================================="
    end subroutine print_matrix

    ! Helper to clean up memory later
    subroutine destroy_matrix(A)
        type(complex_csr_matrix), intent(inout) :: A
        if (allocated(A%row_ptr)) deallocate(A%row_ptr)
        if (allocated(A%col_ind)) deallocate(A%col_ind)
        if (allocated(A%values))  deallocate(A%values)
    end subroutine destroy_matrix


end module complex_sparse_linalg

