module shared_data
    use complex_sparse_linalg
    implicit none
    ! Global parameters
    integer(4) :: ns, Ntot
    real(8), parameter :: epsilon = 1.0d-12 ! epsilon for real comparisons
    real(8) :: dt
    ! Parameter for ADO printouts
    integer(4) :: nprint_ADOs = 100 ! how often to print the ADOs (every nprint_ADOs steps)
    ! other crap
    integer(4) :: Ktot, Imax,  L, lowTcoef_switch
    real(8)::  lowTcoef
    ! Liouvillian
    type(complex_csr_matrix) :: Liouvillian
end module shared_data