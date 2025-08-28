module shared_data
    implicit none
    ! Global parameters
    integer(4) :: Imax, Ktot, ns, L, Ntot
    real(8) :: dt, lowTcoef
    real(8), parameter :: epsilon = 1.0d-12 ! epsilon for real comparisons
    ! Parameter for ADO printouts
    integer(4) :: nprint_ADOs = 100 ! how often to print the ADOs (every nprint_ADOs steps)
    ! Global parameter arrays
    complex(8), allocatable :: gam_ks(:), c_U(:,:), c_D_LEFT(:,:), c_D_RIGHT(:,:)
    integer(4), allocatable :: I0s(:), lengths(:,:), ADO_index(:,:)
    complex(8), allocatable :: s_mat2(:,:), is_mat(:,:), iH_mat(:,:)

end module shared_data