module shared_data
    implicit none
    ! Global parameters
    integer(4) :: ns
    real(8), parameter :: epsilon = 1.0d-12 ! epsilon for real comparisons
    ! Parameter for ADO printouts
    integer(4) :: nprint_ADOs = 100 ! how often to print the ADOs (every nprint_ADOs steps)

end module shared_data