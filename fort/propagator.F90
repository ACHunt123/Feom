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
real(8), parameter  :: Krylov_tol = 1.d-10      ! tolerance for the Krylov subspace
complex(8), public, allocatable :: Krylov_vecs(:,:)     ! Krylov subspace vectors
complex(8) :: L_mat(0:Krylov_dim,0:Krylov_dim)      ! Liouvillian matrix for the Krylov subspace
complex(8) :: ADOs_Krylov(0:Krylov_dim) ! ADOs in the Krylov basis
! Expose some of the subroutines
public :: RK4step, SIAstep, Recalculate_ADOs
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
! Recalculate the ADOs from the Krylov vectors
subroutine Recalculate_ADOs(ADOs)
    implicit none
    complex(8), intent(inout) :: ADOs(Ntot)
    integer(4) :: i,j,k

    !!! reset ADOs
    ADOs=(0.d0,0.d0)
    !!! Recalculate the ADOs from the Krylov basis
    do i=0,Krylov_dim
        do j=1,Ntot
            ADOs(j) = ADOs(j) + ADOs_Krylov(i) * Krylov_vecs(i,j) ! sum over the Krylov vectors
        end do
    end do
    end subroutine Recalculate_ADOs
! Implicitly-Restarted Arnoldi method for the short iterative method
subroutine SIAstep(ADOs)
    implicit none
    complex(8), allocatable :: work(:)
    integer(4) :: info, lwork,i,j
    complex(8) :: L_evals(0:Krylov_dim) ! eigenvalues of the Liouvillian matrix
    complex(8) :: U_Krylov_R(0:Krylov_dim,0:Krylov_dim),U_Krylov_L(0:Krylov_dim,0:Krylov_dim) ! Transformation matrices to diagonalise the Liouvillian in the Krylov subspace
    complex(8), intent(inout) :: ADOs(Ntot)
    complex(8) :: scale, sum
    integer(4) :: k,l

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

    if (abs(ADOs_Krylov(Krylov_dim)).gt.Krylov_tol)    then 
        !!! Recalculate the ADOs from the Krylov basis
        call Recalculate_ADOs(ADOs)

        !!! regenerate the Krylov subspace
        ADOs_Krylov = (0.d0,0.d0)
        ADOs_Krylov(0) = (1.d0,0.d0)


        !!! get the basis vectors for the Krylov subspace
        call generate_krylov_vecs() ! generate the Krylov vectors

        !!! Diagonalise the Liouvillian matrix (DESTROYING IT)
        allocate(work(1)) 
        call zgeev('V','V',Krylov_dim+1,L_mat,Krylov_dim+1,L_evals, &
                U_Krylov_L,Krylov_dim+1,U_Krylov_R,Krylov_dim+1,work,lwork,info)
        if (info.ne.0) stop 'Error in zgeev'
        deallocate(work) 
        ! make the transformation matrices biorthogonal
        do i=0,Krylov_dim
            scale = zdotc(Krylov_dim+1,U_Krylov_L(i,:),1,U_Krylov_R(i,:),1) ! calculate the inner product of the Krylov vectors
            U_Krylov_L(i,:)= U_Krylov_L(i,:)/scale ! scale the left transformation matrix
        end do
        ! Exponentiate eigenvalues with dt
        L_evals(:) = exp(L_evals(:)*dt) ! exponentiate the eigenvalues with the time step
        ! Transform the Liouvillian matrix to the Krylov subspace
        L_mat = 0.d0 ! initialise the Liouvillian matrix in the Krylov subspace to zero
        do k=0,Krylov_dim
            do l=0,Krylov_dim
                sum=(0.d0,0.d0) ! initialise the sum to zero
                do j=0,Krylov_dim
                    sum=sum + U_Krylov_L(k,j) * L_evals(j) * U_Krylov_R(l,j) ! sum over the Krylov vectors
                end do
                L_mat(k,l) = sum ! store the result in the Liouvillian matrix in the Krylov subspace
            end do
        end do
    else
        !!! Propagate the vector with dt
        do k=0,Krylov_dim
            sum=(0.d0,0.d0) ! initialise the sum to zero
            do j=0,Krylov_dim
                sum = sum + L_mat(k,j) * ADOs_Krylov(j) ! sum over the Krylov vectors
            end do
            ADOs_Krylov(k) = sum ! store the result in the ADOs in the Krylov basis
        end do
    end if

    contains
    subroutine generate_krylov_vecs()
        ! This subroutine generates the Krylov vectors for the short iterative method
        implicit none
        integer(4) :: j,k,ni,nj
        real(8) :: norm,beta
        complex(8) :: Phi(Ntot) ! work vector
        

        !Get the zeroth Krylov vector
        L_mat = (0.d0,0.d0)
        norm = dznrm2(Ntot, ADOs, 1)                            ! calculate the norm of the ADOs
        Krylov_vecs(0,:) = ADOs(:)/norm     ! normalise the zeroth Krylov vector
    
        ! calculate the rest of the Krylov vectors
        do k = 0, Krylov_dim-1
            ! calculate the (non-orthonormalised) Krylov vector of the i+1 th order
            ! |phi_{k+1}> =  L|Krylov_vecs(k,:)>
            call get_gradient(Krylov_vecs(k,:),Phi) 

            ! |phi_{k+1}> ==Gram-Schmidt Orthogonalisation==> |Krylov_vecs(k+1,:)>
            do j = 0, k
                ! calculate the inner product of the Krylov vector with the previous Krylov vectors
                beta = zdotc(Ntot, Krylov_vecs(j,:), 1, Phi, 1) ! L_mat(i,j) = <Krylov_vecs(j,:),Phi>
                ! remove the component of the Krylov vector that is in the direction of the previous Krylov vectors
                Phi(:) = Phi(:) - beta*Krylov_vecs(j,:)
                ! store the inner product in the Liouvillian matrix
                L_mat(j,k) = beta   ! L_mat(j,k) = <Krylov_vecs(j,:)|L|Krylov_vecs(k,:)>
            end do
            
            ! normalise the Krylov vector
            norm = dznrm2(Ntot, Phi, 1)                        ! calculate the norm of the Krylov vector
            if (norm.lt.1d-20) then 
                Krylov_vecs(k+1,:) =(0.d0,0.d0) ! if the norm is too small, set the Krylov vector to zero
            else
                Krylov_vecs(k+1,:) = Phi(:)/norm ! normalise the Krylov vector
            end if
        !!!
        end do
        end subroutine generate_krylov_vecs

    end subroutine SIAstep
end module prop_subroutines