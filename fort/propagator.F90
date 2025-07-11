! Fortran module to propagate the HEOM code
module prop_subroutines
use shared_data
use gradient, only: get_gradient
implicit none
! Temporary arrays for RK4 propagation, and constants
complex(8), parameter :: twothirds = 2.d0/3.d0 *(1.d0,0.d0)     ! two thirds, used for the gradient calculation
complex(8), parameter :: third = 1.d0/3.d0* (1.d0,0.d0)         ! third, used for the gradient calculation
complex(8), allocatable :: k1(:,:,:), k2(:,:,:), k3(:,:,:), k4(:,:,:), ktmp(:,:,:),ADOs_tmp(:,:,:), temp_grad(:,:,:)
complex(8) :: dto2,dto6                                         ! half and 1/6 time step, used for the gradient calculation
! Temporary variables for Short Iterative Arnoldi method
integer(4), parameter :: Krylov_dim = 8         ! dimension of the Krylov subspace
real(8), parameter  :: Krylov_tol = 1.d-10      ! tolerance for the Krylov subspace
complex(8), allocatable :: Krylov_vecs(:,:)     ! Krylov subspace vectors
complex(8) :: L_mat(0:Krylov_dim,0:Krylov_dim)      ! Liouvillian matrix for the Krylov subspace
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
    call zscal(Ntot,dt,k3,1) ! scale the gradient by half the time step
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

! Implicitly-Restarted Arnoldi method for the short iterative method
subroutine SAstep(ADOs)
    implicit none
    complex(8), allocatable :: work(:)
    integer(4) :: info, lwork,i,j
    complex(8) :: L_evals(0:Krylov_dim) ! eigenvalues of the Liouvillian matrix
    complex(8) :: L_dt(0:Krylov_dim,0:Krylov_dim) ! Liouvillian matrix in the Krylov subspace for dt
    complex(8) :: U_Krylov_R(0:Krylov_dim,0:Krylov_dim),U_Krylov_L(0:Krylov_dim,0:Krylov_dim) ! Transformation matrices to diagonalise the Liouvillian in the Krylov subspace
    complex(8), intent(inout) :: ADOs(Ntot)

    
    !!! get the basis vectors for the Krylov subspace
    call generate_krylov_vecs() ! generate the Krylov vectors

    !!! Exponentiate the Liouvillian matrix (DESTROYS L_mat)
    ! Diagonalise the Liouvillian matrix (DESTROYING IT)
    allocate(work(1)) 
    call zgeev('V','V',Krylov_dim+1,L_mat,Krylov_dim+1,L_evals, &
               U_Krylov_L,Krylov_dim+1,U_Krylov_R,Krylov_dim+1,work,lwork,info)
    deallocate(work) ! deallocate the work array
    ! Exponentiate eigenvalues with dt
    forall(i=1:Krylov_dim+1) L_evals(i) = exp(-L_evals(i)*dt) ! exponentiate the eigenvalues with the time step
    ! Transform the Liouvillian matrix to the Krylov subspace
    L_dt = 0.d0 ! initialise the Liouvillian matrix in the Krylov subspace to zero
    do i=0:Krylov_dim
        do j=0:Krylov_dim
            L_dt(i,j) = L_evals(i)*U_Krylov_L(i,j) ! L_dt = diag(L_evals) * U_Krylov_L
        end do
    end do
    if (info.ne.0) stop 'Error in zgeev'


    contains
    subroutine generate_krylov_vecs()
        ! This subroutine generates the Krylov vectors for the short iterative method
        implicit none
        integer(4) :: j,k,ni,nj
        real(8) :: norm,beta
        complex(8) :: Phi(Ntot) ! work vector
        interface
            double precision function dznrm2(n, x, incx)
            integer, intent(in) :: n, incx
            complex*16, intent(in) :: x(*)
            end function dznrm2

            complex*16 function zdotc(n, x, incx, y, incy)
            integer, intent(in) :: n, incx, incy
            complex*16, intent(in) :: x(*), y(*)
            end function zdotc
        end interface

        ! initialise the Krylov propagator as zeros
        forall(ni=0:Krylov_dim,nj=0:Krylov_dim) L_mat(ni,nj) = (0.d0,0.d0) ! set the work vector to zero
        
        !get the zeroth Krylov vector
        norm = dznrm2(Ntot, ADOs, 1)                            ! calculate the norm of the ADOs
        forall(ni=1:Ntot) Krylov_vecs(0,ni) = ADOs(ni)/norm     ! normalise the zeroth Krylov vector
    
        ! calculate the rest of the Krylov vectors
        do k = 0, Krylov_dim
            ! calculate the (non-orthonormalised) Krylov vector of the i+1 th order
            call get_gradient(Krylov_vecs(k,:),Phi) 

            do j = 0, k
                ! calculate the inner product of the Krylov vector with the previous Krylov vectors
                beta = zdotc(Ntot, Krylov_vecs(j,:), 1, Phi, 1) ! L_mat(i,j) = <Krylov_vecs(j,:),Phi>
                ! remove the component of the Krylov vector that is in the direction of the previous Krylov vectors
                forall(ni=1:Ntot) Phi(ni) = Phi(ni) - beta*Krylov_vecs(ni,j)
                ! store the inner product in the Liouvillian matrix
                L_mat(k,j) = beta
            end do
            
            ! normalise the Krylov vector
            norm = dznrm2(Ntot, Phi, 1)                        ! calculate the norm of the Krylov vector
            if (norm.lt.1d-20) stop 'Krylov vector norm is too small, stopping iteration'
            forall(ni=1:Ntot) Krylov_vecs(ni,i) = Phi(ni)/norm ! normalise the Krylov vector

        !!!
        end do
        end subroutine generate_krylov_vecs
!
end subroutine SAstep
end module prop_subroutines