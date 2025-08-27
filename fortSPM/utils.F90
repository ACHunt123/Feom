module utils
    implicit none
    interface
        subroutine zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
            implicit none
            character(len=1), intent(in) :: jobu, jobvt
            integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
            integer, intent(out) :: info
            complex*16, intent(inout) :: a(lda, *)
            double precision, intent(out) :: s(*)
            complex*16, intent(out) :: u(ldu, *)
            complex*16, intent(out) :: vt(ldvt, *)
            complex*16, intent(inout) :: work(*)
            double precision, intent(out) :: rwork(*)
        end subroutine zgesvd
    end interface
    contains 
! Inner product function
function innerprod(v1,v2,size)
    ! This function calculates the inner product of two vectors
    ! Input |v1>, |v2> and their size
    ! Output <v1|v2>
    implicit none
    integer(4), intent(in) :: size
    complex(8), intent(in) :: v1(size), v2(size)
    integer(4) :: ii
    complex(8) :: innerprod ! the inner product of the two vectors
    innerprod = (0.d0,0.d0) ! initialise the inner product to zero
    do ii = 1,size
        innerprod = innerprod + conjg(v1(ii)) * v2(ii) ! calculate the inner product
    end do
end function innerprod
! Norm function
function norm(v,size)
    ! This function calculates the norm of a vector
    ! Input |v> and its size
    ! Output sqrt(||v||)
    implicit none
    integer(4), intent(in) :: size
    complex(8), intent(in) :: v(size)
    real(8) :: norm             ! the norm of the vector
    norm = innerprod(v,v,size)  ! calculate the inner product of the vector with itself
    norm = sqrt(real(norm))     ! take the square root of the real
end function norm
! Check the condition number of a matrix
subroutine check_condition_number(A, dim, label)
    implicit none
    integer, intent(in) :: dim
    complex(8), intent(in) :: A(dim, dim)
    character(len=*), intent(in) :: label

    real(8), allocatable :: s(:), rwork(:)
    complex(8), allocatable :: a_copy(:,:), work(:)
    complex(8), allocatable :: u(:,:), vt(:,:)
    integer :: info, lwork

    ! Allocate workspace
    allocate(a_copy(dim, dim))
    a_copy = A

    allocate(s(dim))
    allocate(rwork(5*dim))
    lwork = 10 * dim
    allocate(work(lwork))

    ! Dummy outputs for unused U and VT
    allocate(u(1,1), vt(1,1))

    ! Call LAPACK SVD routine: only singular values needed
    call zgesvd('N', 'N', dim, dim, a_copy, dim, s, u, 1, vt, 1, work, lwork, rwork, info)

    if (info /= 0) then
        print*, 'SVD failed for ', trim(label), ' with info =', info
        stop 'SVD error in check_condition_number'
    end if

    ! Print condition number
    print*, 'Condition number of ', trim(label), ':', maxval(s) / minval(s)

    ! Cleanup
    deallocate(s, rwork, work, a_copy, u, vt)
end subroutine check_condition_number
end module utils