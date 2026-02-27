module gradient
use shared_data, only: Ntot,Liouvillian
use complex_sparse_linalg, only: csr_matvec_complex

implicit none  
! Default everything to private
private
! Expose the gradent function
public :: get_gradient
contains

    ! Function for calculating the gradient of the density, inherits the scope of the vvstep subroutine
subroutine get_gradient(rho,grad)
    complex(8), intent(in) :: rho(Ntot) ! fortran is column major so the last index is the fastest changing
    complex(8), intent(out) :: grad(Ntot)
    call csr_matvec_complex(Liouvillian, rho, grad)

end subroutine get_gradient


end module gradient