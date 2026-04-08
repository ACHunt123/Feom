! Fortran module to propagate the HEOM code
module integrator
use shared_data
use gradient, only: get_gradient
use utils, only: innerprod, norm, check_condition_number
implicit none
! Default everything to private
private
! Temporary arrays for RK4 propagation, and constants
complex(8), parameter :: twothirds = 2.d0/3.d0 *(1.d0,0.d0)     ! two thirds, used for the gradient calculation
complex(8), parameter :: third = 1.d0/3.d0* (1.d0,0.d0)         ! third, used for the gradient calculation
complex(8), allocatable, save :: k1(:), k2(:), k3(:), k4(:), ktmp(:)
complex(8),save :: dto2,dto6                                    ! half and 1/6 time step, used for the gradient calculation
! Temporary variables/Parameters for Short Iterative Arnoldi method
integer(4), parameter :: default_Krylov_dim = 8                 ! default dimension of the Krylov subspace
integer(4), save :: Krylov_dim                                  ! dimension of the Krylov subspace
real(8), parameter  :: Krylov_tol = 1.d-8                       ! tolerance for the Krylov subspace
complex(8), allocatable, save :: Krylov_vecs(:,:)               ! Krylov subspace vectors
complex(8), allocatable, save :: ADOs_Krylov(:),L_mat(:,:)      ! ADOs in the Krylov basis
! Persistent, but private variables for the diagonalization in Arnoldi
complex(8), allocatable, save :: work(:)
real(8), allocatable, save    :: rwork(:)
integer, save                 :: lwork
! Expose some of the subroutines
public :: step_forward, Recalculate_ADOs, init_integrator, cleanup_integrator
! Interface for the matrix diagonalization
interface
        subroutine zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
        character(len=1), intent(in) :: jobvl, jobvr
        integer, intent(in) :: n, lda, ldvl, ldvr, lwork
        complex(8), intent(inout) :: a(lda, *)
        complex(8), intent(out) :: w(*), vl(ldvl, *), vr(ldvr, *), work(*)
        double precision, intent(out) :: rwork(*)
        integer, intent(out) :: info
        end subroutine zgeev
end interface
contains
!!! Initialization of integrators
subroutine init_integrator(ADOs_in)
    complex(8), intent(in) :: ADOs_in(:)
    real(8) :: ADOnorm
    
    Ntot = size(ADOs_in)
#ifdef SIA
    Krylov_dim = min(default_Krylov_dim, Ntot) 
    allocate(Krylov_vecs(Ntot, Krylov_dim)) 
    allocate(ADOs_Krylov(Krylov_dim))
    allocate(L_mat(Krylov_dim, Krylov_dim))

    ! Initialise the subspace state
    ADOnorm = norm(ADOs_in, Ntot)
    if (ADOnorm < 1.d-15) stop 'FATAL: Initial ADOs norm is zero!'
    ADOs_Krylov = (0.0d0, 0.0d0)
    ADOs_Krylov(1) = ADOnorm        ! Keep the physical state intact!
    ADOs_Krylov(Krylov_dim) = 1.0d0 ! Trigger first update
    
    Krylov_vecs = (0.0d0, 0.0d0)
    Krylov_vecs(:, 1) = ADOs_in / ADOnorm 

    ! Call the internal LAPACK workspace setup
    call init_sia_workspace() 
#else
    allocate(k1(Ntot), k2(Ntot), k3(Ntot), k4(Ntot), ktmp(Ntot))
    ! Precomputation of constants
    dto2 = dt / (2.d0, 0.d0)               
    dto6 = dt / (6.d0, 0.d0)               
#endif
end subroutine init_integrator
!!! Cleanup of integrators
subroutine cleanup_integrator()
#ifdef SIA
    if (allocated(Krylov_vecs)) deallocate(Krylov_vecs, ADOs_Krylov, L_mat)
    call cleanup_sia_workspace()
#else
    if (allocated(k1)) deallocate(k1, k2, k3, k4, ktmp)
#endif
end subroutine cleanup_integrator
!!! Forward step (either SIA or RK4)
subroutine step_forward(ADOs)
    implicit none
    complex(8), intent(inout) :: ADOs(Ntot)
#ifdef SIA
            call SIAstep(ADOs)
#else
            call RK4step(ADOs)
#endif
    end subroutine
!!! RK4 subroutines
subroutine RK4step(ADOs)
    implicit none
    complex(8), intent(inout) :: ADOs(Ntot)
    ! Calculate k1
    call get_gradient(ADOs, k1) 
    k1 = k1 * dto2 
    ktmp = ADOs + k1
    ! Calculate k2
    call get_gradient(ktmp, k2) 
    k2 = k2 * dto2
    ktmp = ADOs + k2
    ! Calculate k3
    call get_gradient(ktmp, k3) 
    k3 = k3 * dt
    ktmp = ADOs + k3
    ! Calculate k4
    call get_gradient(ktmp, k4) 
    ! Update the density matrix
    ! Distributing the scalars allows -O3 to fuse the loop -> no temporary arrays 
    ADOs = ADOs + dto6*k4 + twothirds*k2 + third*k3 + third*k1 

    end subroutine RK4step

!!! SIA subroutines
! Recalculate the ADOs from the Krylov vectors
subroutine Recalculate_ADOs(ADOs)
    implicit none
    complex(8), intent(inout) :: ADOs(Ntot)
    ADOs = matmul(Krylov_vecs, ADOs_Krylov)
    end subroutine Recalculate_ADOs
! Initialize the workspace for the diagonalisation of SIA
subroutine init_sia_workspace()
        complex(8) :: query_tmp(1)
        integer :: info
        complex(8) :: dummy_L(Krylov_dim, Krylov_dim)
        complex(8) :: dummy_R(Krylov_dim, Krylov_dim)
        complex(8) :: dummy_evals(Krylov_dim)
        complex(8) :: dummy_mat(Krylov_dim, Krylov_dim)
        ! initialize dummy mat
        dummy_mat = (0.d0, 0.d0)
        ! Pre-allocate rwork (Size is fixed: 2*K)
        if (.not. allocated(rwork)) allocate(rwork(2 * Krylov_dim))
        ! Workspace Query: Ask ZGEEV for the optimal lwork, using dummy arrays to satisfy the interface for the query
        call zgeev('V', 'V', Krylov_dim, dummy_mat, Krylov_dim, dummy_evals, &
                   dummy_L, Krylov_dim, dummy_R, Krylov_dim, &
                   query_tmp, -1, rwork, info)
        lwork = int(real(query_tmp(1)))
        ! Allocate lwork
        if (allocated(work)) deallocate(work)
        allocate(work(lwork))
    end subroutine init_sia_workspace
! Implicitly-Restarted Arnoldi method for the short iterative method
subroutine SIAstep(ADOs)
    implicit none
    integer :: info,i,j,k
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
        call zgeev('V','V', Krylov_dim, L_mat, Krylov_dim, L_evals, &
                U_Krylov_L, Krylov_dim, U_Krylov_R, Krylov_dim, work, lwork, rwork, info)
        if (info /= 0) stop 'Error in zgeev'
        ! make the transformation matrices biorthogonal
        do i=1,Krylov_dim
            scale = innerprod(U_Krylov_L(:,i) ,U_Krylov_R(:,i),Krylov_dim) ! calculate the inner product of the Krylov vectors
            U_Krylov_L(:,i)= U_Krylov_L(:,i)/dconjg(scale) ! scale the left transformation matrix
        end do

        !!! Exponentiate eigenvalues with dt to get the eigenvalues of the propagator
        L_evals = zexp(L_evals*dt) ! exponentiate the eigenvalues with the time step

        if (.false.) then !Tests: BIORTHOGONALITY, CONDITION NUMBER, EIGENVALUES
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
        integer(4) :: j,k
        real(8) :: ADOnorm
        complex(8) :: Phi(Ntot),beta ! work vector, projection
        
        !Get the zeroth Krylov vector
        L_mat = (0.d0,0.d0)
        ADOnorm = norm(ADOs,Ntot)                            ! calculate the norm of the ADOs
        Krylov_vecs= (0.d0,0.d0) ! initialise the Krylov vectors to zero
        phi = (0.d0,0.d0) ! initialise the work vector to zero
        Krylov_vecs(:,1) = ADOs(:)/ADOnorm     ! normalise the zeroth Krylov vector
    
        ! calculate the rest of the Krylov vectors
        do k = 1, Krylov_dim
            ! calculate the (non-orthonormalised) Krylov vector of the i+1 th order
            ! |phi_{k+1}> =  L|Krylov_vecs(:,k)>
            call get_gradient(Krylov_vecs(:,k),Phi) 

            ! |phi_{k+1}> ==Gram-Schmidt Orthogonalisation==> |Krylov_vecs(:,k+1)>
            do j = 1, k
                ! calculate the inner product of the Krylov vector with the previous Krylov vectors
                beta = innerprod(Krylov_vecs(:,j),Phi, Ntot) ! L_mat(i,j) = <Krylov_vecs(:,j),Phi>
                ! remove the component of the Krylov vector that is in the direction of the previous Krylov vectors
                Phi = Phi - beta*Krylov_vecs(:,j)
                ! store the inner product in the Liouvillian matrix
                L_mat(j,k) = beta   ! L_mat(j,k) = <Krylov_vecs(:,j)|L|Krylov_vecs(:,k)>
                
            end do
            
            ! normalise the Krylov vector, and store the k+1,k th element of the Liouvillian matrix
            ADOnorm = norm(Phi,Ntot)                        ! calculate the norm of the Krylov vector
            if (k.lt.Krylov_dim) then
                L_mat(k+1,k) = ADOnorm ! <Krylov_vecs(:,k+1)|L|Krylov_vecs(:,k)> = <phi| * (|phi> + stuff orthogal to |phi>) = <phi|phi>
                if (abs(ADOnorm).lt.1d-20) then 
                    Krylov_vecs(:,k+1) =(0.d0,0.d0) ! if the norm is too small, set the Krylov vector to zero
                else
                    Krylov_vecs(:,k+1) = Phi/ADOnorm ! normalise the Krylov vector
                end if
            endif
        !!!
        end do
        end subroutine generate_krylov_vecs

    end subroutine SIAstep
! Initialize the workspace for the diagonalisation of SIA
subroutine cleanup_sia_workspace()
    ! Check if allocated before deallocating to avoid crashes
    if (allocated(work)) then
        deallocate(work)
    endif
    if (allocated(rwork)) then
        deallocate(rwork)
    endif
    lwork = 0
    end subroutine cleanup_sia_workspace
end module integrator