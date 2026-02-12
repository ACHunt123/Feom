module shared_data
    implicit none
    ! Global parameters
    integer(4) :: Imax, Ktot, ns, L, Ntot, hbar
    real(8) :: dt
    real(8), parameter :: epsilon = 1.0d-12 ! epsilon for real comparisons
    ! Global parameter arrays
    complex(8), allocatable :: gam_ks(:), c_U(:,:), c_D_LEFT(:,:), c_D_RIGHT(:,:)
    integer(4), allocatable :: I0s(:), lengths(:,:), ADO_index(:,:)
    complex(8), allocatable :: is_mat(:,:), iH_mat(:,:)
    ! terminator superoperator
    complex(8), allocatable :: Xi_terminator(:,:,:)

end module shared_data