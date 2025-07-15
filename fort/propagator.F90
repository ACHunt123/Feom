! Fortran module to propagate the HEOM code
module prop_subroutines
use shared_data
use gradient, only: get_gradient
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
real(8), parameter  :: Krylov_tol = 1.d-10              ! tolerance for the Krylov subspace
complex(8) :: L_mat(Krylov_dim,Krylov_dim)          ! Liouvillian matrix for the Krylov subspace
complex(8), public, allocatable :: Krylov_vecs(:,:)     ! Krylov subspace vectors
complex(8), public :: ADOs_Krylov(Krylov_dim)         ! ADOs in the Krylov basis
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


!    ! Function for killing off ADOs that have a small norm
subroutine changeN(ADOs,active,active0,Nactive,Nactive0)
    implicit none
    complex(8), intent(inout), allocatable :: ADOs(:,:,:)
    integer(4), intent(inout) :: active(Imax), active0(Imax), Nactive, Nactive0
    integer(4) :: I, si, sj
    real(8) :: norm
    
    ! reallocate the temporary array if required
    if (Nactive0.ne.Nactive) then
        deallocate(ADOs_tmp)
        allocate(ADOs_tmp(Nactive,ns,ns))
        ADOs_tmp = 0.d0
    end if

    ! reset number of active ADOs (as we are about to kill the small ones)
    Nactive = 0
    do I = 1, Imax
        ! all ADOS
        if (active(I).eq.0) cycle !skip if ADO not active at the end of the timestep
        ! all ADOs that are active at end of timestep
        if (active0(I).eq.0) then !if the ADO was not active at the start of the timestep, then it is not active now
            Nactive = Nactive + 1
            active(I) = Nactive
            ADOs_tmp(Nactive,:,:) =  0.d0 ! create new ADO in list
            cycle
        endif
        ! all ados that were active at start and end of timestep (now we kill off some)
        norm = 0.d0
        do si = 1, ns
            do sj = 1, ns
                norm = norm + abs(ADOs(active0(I),si,sj))**2
            end do
        end do
        if (norm.gt.tolerance) then
            Nactive = Nactive + 1
            active(I) = Nactive
            ADOs_tmp(Nactive,:,:) = ADOs(active0(I),:,:)
            cycle
        end if
        ! the only left are the ones that are too small
        active(I) = 0
        
    end do

    if (Nactive0.ne.Nactive) then
        deallocate(ADOs)
        allocate(ADOs(Nactive,ns,ns))
    end if
    do I = 1, Nactive
        ADOs(I,:,:) = ADOs_tmp(I,:,:)
    end do
    
    end subroutine changeN


! rk4 propagation of HEOM for one step
subroutine RK4step(ADOs)
    implicit none
    complex(8), allocatable, intent(inout) :: ADOs(:,:,:)
    if(Imax.gt.2147483647) stop 'Imax is too large for the ADO index array'
    ! Precomputation of constants
    dto2 = dt/(2.d0,0.d0)               ! half the time step, used for the gradient calculation
    dto6 = dt/(6.d0,0.d0)               ! 1/6 the time step, used for the gradient calculation

    ! reallovate the temporary arrays if required
    #ifdef Prune
    if (Nactive0.ne.Nactive) then
        deallocate(k1,k2,k3,k4,ktmp,temp_grad)
        allocate(k1(Nactive,ns,ns), k2(Nactive,ns,ns), k3(Nactive,ns,ns), k4(Nactive,ns,ns), ktmp(Nactive,ns,ns),temp_grad(Nactive,ns,ns))
    end if
    #endif

    ! initialise nactive0 and active0
    ! those without 0 will be changed during propagations, those with 0 will not
    Nactive0 = Nactive
    active0 = active

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

    
    ! Kill off ADOs with small norm
    #ifdef Prune
    call changeN(ADOs,active,active0,Nactive,Nactive0)
    #endif
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

    if (abs(ADOs_Krylov(Krylov_dim)).gt.Krylov_tol)    then 
        ! print*, 'Recalculating Krylov subspace...'
        !!! Recalculate the ADOs from the Krylov basis
        call Recalculate_ADOs(ADOs)

        !!! regenerate the Krylov subspace using these new ADOs
        ADOs_Krylov = (0.d0,0.d0)
        ADOs_Krylov(1) = (1.d0,0.d0)

        ! print*, 'Generating Krylov subspace...'
        !!! get the basis vectors for the Krylov subspace
        call generate_krylov_vecs(ADOs) ! generate the Krylov vectors
        L_mat_copy = L_mat ! make a copy of the Liouvillian matrix
        !!! Diagonalise the Liouvillian matrix (DESTROYING IT)
        lwork= 2*(Krylov_dim)*(Krylov_dim) ! workspace size for the LAPACK routine
        allocate(work(lwork)) ! allocate the work array for the LAPACK routine
        allocate(rwork(2*(Krylov_dim))) ! allocate the real work array for the LAPACK routine
        call zgeev('V','V', Krylov_dim, L_mat, Krylov_dim, L_evals, &
                U_Krylov_L, Krylov_dim, U_Krylov_R, Krylov_dim, work, lwork, rwork, info)
        deallocate(work) ! deallocate the work array
        deallocate(rwork) ! deallocate the real work array
        ! print*, 'Eigenvalues of the Liouvillian matrix:'
        ! do i=1,Krylov_dim
        !     print*, real(L_evals(i)), aimag(L_evals(i))
        ! end do
        if (info /= 0) stop 'Error in zgeev'
        ! Test that the eigenvectors are actually eigenvectors
        do i = 1,Krylov_dim
            vec_R=U_Krylov_R(:,i)
            vec_L=conjg(U_Krylov_L(:,i))
            lambda=L_evals(i)
            Av= matmul(L_mat_copy,vec_R)
            vA= matmul(vec_L,L_mat_copy)
            ! print*,sum(abs(Av-lambda*vec_R))
            ! print*,sum(abs(vA-lambda*vec_L))
        end do
        ! print*, 'Krylov subspace generated.'

        ! make the transformation matrices biorthogonal
        ! print*, 'Biorthogonalising the transformation matrices...'
        do i=1,Krylov_dim
            scale = innerprod(U_Krylov_L(:,i) ,U_Krylov_R(:,i),Krylov_dim) ! calculate the inner product of the Krylov vectors
            U_Krylov_L(:,i)= U_Krylov_L(:,i)/dconjg(scale) ! scale the left transformation matrix
        end do
        ! test that the transformation matrices are biorthogonal
        do i=1,Krylov_dim
            do j=1,Krylov_dim
                ! print*,abs(innerprod(U_Krylov_L(:,i),U_Krylov_R(:,j),Krylov_dim))
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

        ! print*, 'exponentiating the eigenvalues...'
        ! Exponentiate eigenvalues with dt
        L_evals = zexp(L_evals*dt) ! exponentiate the eigenvalues with the time step

        ! Transform the Liouvillian matrix to the Krylov subspace
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
       
        ! print*, 'Liouvillian matrix in the Krylov subspace:'
        ! do i=1,Krylov_dim
        !     write(*,'(10F10.5)') real(L_mat(i,:))
        ! end do
        ! print*, 'copy of the original'
        ! do i=1,Krylov_dim
        !     write(*,'(10F10.5)') real(L_mat_copy(i,:))
        ! end do

        
        ! stop
    else
        !!! Propagate the vector with dt
        do k=1,Krylov_dim
            csum=(0.d0,0.d0) ! initialise the csum to zero
            do j=1,Krylov_dim
                csum = csum + L_mat(k,j) * ADOs_Krylov(j) ! sum over the Krylov vectors
            end do
            ADOs_Krylov(k) = csum ! store the result in the ADOs in the Krylov basis
        end do
    end if

    contains
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

    subroutine generate_krylov_vecs(ADOs) 
        ! This subroutine generates the Krylov vectors for the short iterative method
        implicit none
        complex(8), intent(in) :: ADOs(Ntot) ! the ADOs in the Krylov basis
        integer(4) :: j,k,ni,nj
        real(8) :: norm,beta
        complex(8) :: Phi(Nactive0,ns,ns) ! work vector
        
        !Get the zeroth Krylov vector
        L_mat = (0.d0,0.d0)
        norm = dznrm2(Ntot, ADOs, 1)                            ! calculate the norm of the ADOs
        Krylov_vecs= (0.d0,0.d0) ! initialise the Krylov vectors to zero
        phi = (0.d0,0.d0) ! initialise the work vector to zero
        Krylov_vecs(1,:) = ADOs(:)/norm     ! normalise the zeroth Krylov vector
    
        ! calculate the rest of the Krylov vectors
        do k = 1, Krylov_dim
            ! calculate the (non-orthonormalised) Krylov vector of the i+1 th order
            ! |phi_{k+1}> =  L|Krylov_vecs(k,:)>
            call get_gradient(reshape(Krylov_vecs(k,:),[Nactive0,ns,ns]),Phi) 

            ! |phi_{k+1}> ==Gram-Schmidt Orthogonalisation==> |Krylov_vecs(k+1,:)>
            do j = 1, k
                ! calculate the inner product of the Krylov vector with the previous Krylov vectors
                beta = innerprod(Krylov_vecs(j,:),reshape(Phi,[Ntot]), Ntot) ! L_mat(i,j) = <Krylov_vecs(j,:),Phi>
                ! remove the component of the Krylov vector that is in the direction of the previous Krylov vectors
                Phi = Phi - beta*reshape(Krylov_vecs(j,:),[Nactive0,ns,ns])
                ! store the inner product in the Liouvillian matrix
                L_mat(j,k) = beta   ! L_mat(j,k) = <Krylov_vecs(j,:)|L|Krylov_vecs(k,:)>
                
            end do
            
            ! normalise the Krylov vector, and store the k+1,k th element of the Liouvillian matrix
            norm = dznrm2(Ntot, Phi, 1)                        ! calculate the norm of the Krylov vector
            if (k.lt.Krylov_dim) then
                L_mat(k+1,k) = norm ! <Krylov_vecs(k+1,:)|L|Krylov_vecs(k,:)> = <phi| * (|phi> + stuff orthogal to |phi>) = <phi|phi>
                if (abs(norm).lt.1d-20) then 
                    Krylov_vecs(k+1,:) =(0.d0,0.d0) ! if the norm is too small, set the Krylov vector to zero
                else
                    Krylov_vecs(k+1,:) = reshape(Phi,[Ntot])/norm ! normalise the Krylov vector
                end if
            endif
        !!!
        end do
        end subroutine generate_krylov_vecs

    end subroutine SIAstep
end module prop_subroutines