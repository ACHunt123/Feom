module shared_data
    implicit none
    ! Global parameters
    integer(4) :: lowTcoef_switch
    integer(4) :: Imax, Ktot, ns, L
    real(8) :: dt, lowTcoef
    real(8), parameter :: epsilon = 1.0d-12 ! epsilon for real comparisons
    real(8), parameter :: tolerance = 1.d-20! tolerance for killing off ADOs with small norm
    ! Parameter for ADO printouts
    integer(4) :: nprint_ADOs = 100 ! how often to print the ADOs (every nprint_ADOs steps)
    ! Global arrays used for allocation reduction
    integer(4), allocatable :: active(:), active0(:) ! active[0] is to store list of active ADOs [at start of timestep] and their indices
    integer(4) :: Nactive0, Nactive ! Each vvstep, arrays are allocated (Nactive0,ns,ns), with the index of rho(I,:,:) being rho(active0(I),:,:)

    ! Global parameter arrays
    complex(8), allocatable :: gam_ks(:), c_U(:,:), c_D_LEFT(:,:), c_D_RIGHT(:,:)
    integer(4), allocatable :: I0s(:), lengths(:,:), ADO_index(:,:)
    complex(8), allocatable :: s_mat2(:,:), is_mat(:,:), iH_mat(:,:)

end module shared_data