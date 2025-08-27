! Fortran module to propagate the HEOM code
module prop_subroutines
use shared_data
use gradient, only: get_gradient
use utils, only: innerprod,norm, check_condition_number
implicit none
! Default everything to private
private
! Temporary arrays for RK4 propagation, and constants
complex(8), parameter :: twothirds = 2.d0/3.d0 *(1.d0,0.d0)     ! two thirds, used for the gradient calculation
complex(8), parameter :: third = 1.d0/3.d0* (1.d0,0.d0)         ! third, used for the gradient calculation
complex(8), public, allocatable :: k1(:,:,:), k2(:,:,:), k3(:,:,:), k4(:,:,:), ktmp(:,:,:),ADOs_tmp(:,:,:), temp_grad(:,:,:)
complex(8) :: dto2,dto6                                         ! half and 1/6 time step, used for the gradient calculation
! Temporary variables for Short Iterative Arnoldi method
integer(4), public, parameter :: Krylov_dim = 8         ! dimension of the Krylov subspace
real(8), parameter  :: Krylov_tol = 1.d-8               ! tolerance for the Krylov subspace
complex(8) :: L_mat(Krylov_dim,Krylov_dim)              ! Liouvillian matrix for the Krylov subspace
complex(8), public, allocatable :: Krylov_vecs(:,:)     ! Krylov subspace vectors
complex(8), public :: ADOs_Krylov(Krylov_dim)           ! ADOs in the Krylov basis
! Expose some of the subroutines
public :: RK4step, SIAstep, Recalculate_ADOs
! Interfaces for the LAPACK routines
interface
        double precision function dznrm2(n, x, incx)
        integer, intent(in) :: n, incx
        complex*16, intent(in) :: x(*)
        end function dznrm2

        complex*16 function zdotc(n, x, incx, y, incy)
        integer, intent(in) :: n, incx, incy
        complex*16, intent(in) :: x(*), y(*)
        end function zdotc

        subroutine zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
        character(len=1), intent(in) :: jobvl, jobvr
        integer, intent(in) :: n, lda, ldvl, ldvr, lwork
        complex(8), intent(inout) :: a(lda, *)
        complex(8), intent(out) :: w(*), vl(ldvl, *), vr(ldvr, *), work(*)
        double precision, intent(out) :: rwork(*)
        integer, intent(out) :: info
        end subroutine zgeev

        subroutine zcopy(n, x, incx, y, incy)
        integer, intent(in) :: n, incx, incy
        complex(8), intent(in) :: x(*)
        complex(8), intent(out) :: y(*)
        end subroutine zcopy

        subroutine zaxpy(n, alpha, x, incx, y, incy)
        integer, intent(in) :: n, incx, incy
        complex(8), intent(in) :: alpha, x(*)
        complex(8), intent(inout) :: y(*)
        end subroutine zaxpy

        subroutine zscal(n, alpha, x, incx)
        integer, intent(in) :: n, incx
        complex(8), intent(in) :: alpha
        complex(8), intent(inout) :: x(*)
        end subroutine zscal
end interface
contains


! rk4 propagation of HEOM for one step
subroutine RK4step(ADOs)
    implicit none
    complex(8), intent(inout) :: ADOs(Imax,ns,ns)
    if(Imax.gt.2147483647) stop 'Imax is too large for the ADO index array'
    ! Precomputation of constants
    dto2 = dt/(2.d0,0.d0)               ! half the time step, used for the gradient calculation
    dto6 = dt/(6.d0,0.d0)               ! 1/6 the time step, used for the gradient calculation

    ! Calculate the k values for the Runge-Kutta method
    call get_gradient(ADOs,k1) ! calculate the gradient of the density matrix
    ! k1 = temp_grad*dto2 !need to split the computation here into two lines to avoid a bug
    call zscal(Ntot,dto2,k1,1) ! scale the gradient by half the time step
    ! ktmp = ADOs+k1
    call zcopy(Ntot, ADOs, 1, ktmp, 1)              ! copy the ADOs to a temporary array
    call zaxpy(Ntot, (1.d0,0.d0), k1, 1, ktmp, 1)   ! add the scaled gradient to the ADOs

    call get_gradient(ktmp,k2) ! calculate the gradient of the updated density matrix
    ! k2 = temp_grad*dto2
    call zscal(Ntot,dto2,k2,1) ! scale the gradient by half the time step
    ! ktmp = ADOs+k2
    call zcopy(Ntot, ADOs, 1, ktmp, 1)              ! copy the ADOs to a temporary array
    call zaxpy(Ntot, (1.d0,0.d0), k2, 1, ktmp, 1)   ! add the scaled gradient to the ADOs

    call get_gradient(ktmp,k3) ! calculate the gradient of the updated density matrix
    ! k3 = temp_grad*dt
    call zscal(Ntot,dt*(1.d0,0.d0),k3,1) ! scale the gradient by half the time step
    ! ktmp = ADOs+k3
    call zcopy(Ntot, ADOs, 1, ktmp, 1)              ! copy the ADOs to a temporary array
    call zaxpy(Ntot, (1.d0,0.d0), k3, 1, ktmp, 1)   ! add the scaled gradient to the ADOs

    call get_gradient(ktmp,k4) 

    ! ! Update the density matrix
    ! ADOs = ADOs + dto6*k4 + twothirds*k2 + (k3 + k1)*third 
    call zaxpy(Ntot, dto6, k4, 1, ADOs, 1)      ! ADOs = ADOs + dto6 * k4
    call zaxpy(Ntot, twothirds, k2, 1, ADOs, 1) ! ADOs = ADOs + twothirds * k2
    call zaxpy(Ntot, third, k3, 1, ADOs, 1)     ! ADOs = ADOs + third * k3
    call zaxpy(Ntot, third, k1, 1, ADOs, 1)     ! ADOs = ADOs + third * k1

    end subroutine
! Recalculate the ADOs from the Krylov vectors
subroutine Recalculate_ADOs(ADOs)
    implicit none
    complex(8), intent(inout) :: ADOs(Ntot)
    integer(4) :: i,j
    !!! reset ADOs
    ADOs=(0.d0,0.d0)
    !!! Recalculate the ADOs from the Krylov basis
    do i=1,Krylov_dim
        do j=1,Ntot
            ADOs(j) = ADOs(j) + ADOs_Krylov(i) * Krylov_vecs(i,j) ! sum over the Krylov vectors
        end do
    end do
    end subroutine Recalculate_ADOs
! Implicitly-Restarted Arnoldi method for the short iterative method
subroutine SIAstep(ADOs)
    implicit none
    complex(8), allocatable :: work(:)
    double precision, allocatable :: rwork(:)
    integer :: info, lwork,i,j,k
    complex(8) :: L_evals(Krylov_dim) ! eigenvalues of the Liouvillian matrix
    complex(8) :: U_Krylov_R(Krylov_dim,Krylov_dim),U_Krylov_L(Krylov_dim,Krylov_dim) ! Transformation matrices to diagonalise the Liouvillian in the Krylov subspace
    complex(8), intent(inout) :: ADOs(Ntot)
    complex(8) :: scale, csum ! temporary variables
    complex(8) :: vec_R(Krylov_dim),vec_L(Krylov_dim),lambda,Av(Krylov_dim),vA(Krylov_dim),L_mat_copy(Krylov_dim,Krylov_dim) ! temporary variables for the eigenvector and eigenvalue

    ! get the basis vectors for the Krylov subspace (if it is needed)
    if (abs(ADOs_Krylov(Krylov_dim)).gt.Krylov_tol*norm(ADOs_Krylov,Krylov_dim))    then 

        !!! Recalculate the ADOs from the Krylov basis
        call Recalculate_ADOs(ADOs)

        !!! regenerate the Krylov subspace using these new ADOs
        
        !!! get the basis vectors for the Krylov subspace
        call generate_krylov_vecs(ADOs)         ! generate the Krylov vectors
        ADOs_Krylov = (0.d0,0.d0)               ! in this new basis, the first vector is the ADOs/norm
        ADOs_Krylov(1) = norm(ADOs,Ntot)  ! so the coefficient is the norm of the ADOs
        L_mat_copy = L_mat ! make a copy of the Liouvillian matrix

        !!! Diagonalise the Liouvillian matrix (DESTROYING IT)
        lwork= 2*(Krylov_dim)*(Krylov_dim) ! workspace size for the LAPACK routine
        allocate(work(lwork)) ! allocate the work array for the LAPACK routine
        allocate(rwork(2*(Krylov_dim))) ! allocate the real work array for the LAPACK routine
        call zgeev('V','V', Krylov_dim, L_mat, Krylov_dim, L_evals, &
                U_Krylov_L, Krylov_dim, U_Krylov_R, Krylov_dim, work, lwork, rwork, info)
        deallocate(work) ! deallocate the work array
        deallocate(rwork) ! deallocate the real work array
        if (info /= 0) stop 'Error in zgeev'
        ! make the transformation matrices biorthogonal
        do i=1,Krylov_dim
            scale = innerprod(U_Krylov_L(:,i) ,U_Krylov_R(:,i),Krylov_dim) ! calculate the inner product of the Krylov vectors
            U_Krylov_L(:,i)= U_Krylov_L(:,i)/dconjg(scale) ! scale the left transformation matrix
        end do

        !!! Exponentiate eigenvalues with dt to get the eigenvalues of the propagator
        L_evals = zexp(L_evals*dt) ! exponentiate the eigenvalues with the time step

        if (0) then !Tests: BIORTHOGONALITY, CONDITION NUMBER, EIGENVALUES
            ! eigenvalues and eigenvectors are correct
            do i = 1,Krylov_dim
                vec_R=U_Krylov_R(:,i)
                vec_L=conjg(U_Krylov_L(:,i))
                lambda=L_evals(i)
                Av= matmul(L_mat_copy,vec_R)
                vA= matmul(vec_L,L_mat_copy)
                if (sum(abs(Av-lambda*vec_R))>1e-10) print*, 'Error in eigenvector calculation for right vector i=',i, sum(abs(Av-lambda*vec_R))
                if (sum(abs(vA-lambda*vec_L))>1e-10) print*, 'Error in eigenvector calculation for left vector i=',i, sum(abs(vA-lambda*vec_L))
            end do
            ! test that the transformation matrices are biorthogonal
            do i=1,Krylov_dim
                do j=1,Krylov_dim
                    if (i.eq.j) then
                        if ((abs(innerprod(U_Krylov_L(:,i),U_Krylov_R(:,j),Krylov_dim))-1.d0) > 1.d-10) then
                            print*, 'Biorthogonalisation failed for i=j=',i,abs(innerprod(U_Krylov_L(i,:),U_Krylov_R(j,:),Krylov_dim))
                            stop
                        end if
                    else
                        if (abs(innerprod(U_Krylov_L(:,i),U_Krylov_R(:,j),Krylov_dim)) > 1.d-10) then
                            print*, 'Biorthogonalisation failed for i=',i,' j=',j,abs(innerprod(U_Krylov_L(i,:),U_Krylov_R(j,:),Krylov_dim))
                            stop
                        end if
                    end if
                end do
            end do
            ! test how conditioned the eigenvalues are
            call check_condition_number(U_Krylov_R, Krylov_dim, 'U_Krylov_R')
            call check_condition_number(U_Krylov_L, Krylov_dim, 'U_Krylov_L')
            ! test that the eigenvalues are not too big
            if (any(abs(L_evals).gt.100.d0)) stop 'Error: eigenvalues of the Liouvillian matrix are too large'
        end if

        !!! Transform the Liouvillian matrix to the Krylov subspace
        L_mat = 0.d0 ! initialise the Liouvillian matrix in the Krylov subspace to zero
        do i=1,Krylov_dim
            do j=1,Krylov_dim
                csum=(0.d0,0.d0) ! initialise the csum to zero
                do k=1,Krylov_dim
                    csum=csum + U_Krylov_R(i,k)* L_evals(k) * dconjg(U_Krylov_L(j,k))  ! sum over the Krylov vectors
                end do
                L_mat(i,j) = csum ! store the result in the Liouvillian matrix in the Krylov subspace
            end do
        end do
       
    end if

    !!! Propagate the vector with dt 
    ADOs_Krylov = matmul(L_mat,ADOs_Krylov) ! propagate the ADOs in the Krylov basis with the Liouvillian matrix
    

    contains
    subroutine generate_krylov_vecs(ADOs) 
        ! This subroutine generates the Krylov vectors for the short iterative method
        implicit none
        complex(8), intent(in) :: ADOs(Ntot) ! the ADOs in the Krylov basis
        integer(4) :: j,k,ni,nj
        real(8) :: ADOnorm,beta
        complex(8) :: Phi(Imax,ns,ns) ! work vector
        
        !Get the zeroth Krylov vector
        L_mat = (0.d0,0.d0)
        ADOnorm = norm(ADOs,Ntot)                            ! calculate the norm of the ADOs
        Krylov_vecs= (0.d0,0.d0) ! initialise the Krylov vectors to zero
        phi = (0.d0,0.d0) ! initialise the work vector to zero
        Krylov_vecs(1,:) = ADOs(:)/ADOnorm     ! normalise the zeroth Krylov vector
    
        ! calculate the rest of the Krylov vectors
        do k = 1, Krylov_dim
            ! calculate the (non-orthonormalised) Krylov vector of the i+1 th order
            ! |phi_{k+1}> =  L|Krylov_vecs(k,:)>
            call get_gradient(reshape(Krylov_vecs(k,:),[Imax,ns,ns]),Phi) 

            ! |phi_{k+1}> ==Gram-Schmidt Orthogonalisation==> |Krylov_vecs(k+1,:)>
            do j = 1, k
                ! calculate the inner product of the Krylov vector with the previous Krylov vectors
                beta = innerprod(Krylov_vecs(j,:),reshape(Phi,[Ntot]), Ntot) ! L_mat(i,j) = <Krylov_vecs(j,:),Phi>
                ! remove the component of the Krylov vector that is in the direction of the previous Krylov vectors
                Phi = Phi - beta*reshape(Krylov_vecs(j,:),[Imax,ns,ns])
                ! store the inner product in the Liouvillian matrix
                L_mat(j,k) = beta   ! L_mat(j,k) = <Krylov_vecs(j,:)|L|Krylov_vecs(k,:)>
                
            end do
            
            ! normalise the Krylov vector, and store the k+1,k th element of the Liouvillian matrix
            ADOnorm = norm(Phi,Ntot)                        ! calculate the norm of the Krylov vector
            if (k.lt.Krylov_dim) then
                L_mat(k+1,k) = ADOnorm ! <Krylov_vecs(k+1,:)|L|Krylov_vecs(k,:)> = <phi| * (|phi> + stuff orthogal to |phi>) = <phi|phi>
                if (abs(ADOnorm).lt.1d-20) then 
                    Krylov_vecs(k+1,:) =(0.d0,0.d0) ! if the norm is too small, set the Krylov vector to zero
                else
                    Krylov_vecs(k+1,:) = reshape(Phi,[Ntot])/ADOnorm ! normalise the Krylov vector
                end if
            endif
        !!!
        end do
        end subroutine generate_krylov_vecs

    end subroutine SIAstep
end module prop_subroutines