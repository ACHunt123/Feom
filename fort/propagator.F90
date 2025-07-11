! Fortran module to propagate the HEOM code
module prop_subroutines
use shared_data, only: Imax, ns, dt, tolerance, Nactive0, Nactive, active, active0
use gradient
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
complex(8) :: L_mat(Krylov_dim,Krylov_dim)      ! Liouvillian matrix for the Krylov subspace
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
    complex(8), intent(in) :: ADOs(Ntot,ns,ns)
    
    
end subroutine SAstep
end module prop_subroutines