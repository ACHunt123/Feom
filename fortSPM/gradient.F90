module gradient
use shared_data, only: Ntot

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

    ! ! Loop over the ADOs
    ! do I = 1, Imax
    !     rhoI = rho(I,:,:) ! temporary variable for the ADO (same done for the gradient, but needn't be initiallised as is overwritten)

    !     ! Get the n values for the ADOs
    !     n_ks = ADO_index(I,:)

    !     ! Execute the diagonal superoperator terms
    !     #ifdef USEZGEMM
    !     call ZGEMM('N','N',ns,ns,ns,-(1.0d0,0.0d0),iH_mat,ns,rhoI,ns,(0.0d0,0.0d0),gradI,ns) ! gradI = -iH_mat * rhoI
    !     call ZGEMM('N','N',ns,ns,ns,(1.0d0,0.0d0),rhoI,ns,iH_mat,ns,(1.0d0,0.0d0),gradI,ns)  ! gradI = gradI + rhoI * iH_mat
    !     gradI = gradI  - sum(n_ks * gam_ks) * rhoI
    !     #else
    !     gradI = ( - matmul(iH_mat,rhoI) + matmul(rhoI,iH_mat)) &
    !             - sum(n_ks * gam_ks) * rhoI 
    !     #endif

    !     ! Itziki Trucation (if present)
    !     #ifdef LowTCorr
    !         gradI = gradI + lowTcoef  &
    !         * ( matmul(s_mat2,rhoI) + matmul(rhoI,s_mat2) + 2.d0*matmul(matmul(is_mat,rhoI),is_mat)) !note the + on last term is as (is)Rho(is) = - 2 s Rho s
    !     #endif

    !     ! Execute the off-diagonal superoperator terms
    !     !NOTE: The [commented lines] are what the code is really doing (we put in the precalculated values to speed up the code)
    !     do ki = 1, (Ktot)
    !         nk = n_ks(ki)
    !         ! Calculate the indices of the ADOs with nk+1 and nk-1
    !         I_nkp1 = I_nk_plusminus(I,ki,+1)
    !         I_nkm1 = I_nk_plusminus(I,ki,-1)

    !         if (I_nkp1.ne.-1) then !check if the index is valid

    !             rhoInkp1 = rho(I_nkp1,:,:)
    !             #ifdef USEZGEMM
    !             call ZGEMM('N','N',ns,ns,ns,c_U(ki,nk),rhoInkp1,ns,is_mat,ns,(1.0d0,0.0d0),gradI,ns)    ! gradI = gradI + c_U(ki,nk) * rhoInkp1 * is_mat
    !             call ZGEMM('N','N',ns,ns,ns,-c_U(ki,nk),is_mat,ns,rhoInkp1,ns,(1.0d0,0.0d0),gradI,ns)   ! gradI = gradI - c_U(ki,nk) * is_mat * rhoInkp1
    !             #else
    !             gradI = gradI &
    !             +  c_U(ki,nk) * ( - matmul(is_mat,rhoInkp1) + matmul(rhoInkp1,is_mat))
    !             !   + [sqrt((nk+1)*abs(Ck)) * ( - matmul(is_mat,rho(I_nkp1,:,:)) + matmul(rho(I_nkp1,:,:),is_mat))]
    !             #endif
    !         endif


    !         if (I_nkm1.ne.-1) then
               
    !             rhoInkm1 = rho(I_nkm1,:,:)
    !             #ifdef USEZGEMM
    !             call ZGEMM('N','N',ns,ns,ns,c_D_LEFT(ki,nk),is_mat,ns,rhoInkm1,ns,(1.0d0,0.0d0),gradI,ns)  ! gradI = gradI + c_D_LEFT(ki,nk) * is_mat * rhoInkm1
    !             call ZGEMM('N','N',ns,ns,ns,c_D_RIGHT(ki,nk),rhoInkm1,ns,is_mat,ns,(1.0d0,0.0d0),gradI,ns) ! gradI = gradI + c_D_RIGHT(ki,nk) * rhoInkm1 * is_mat
    !             #else
    !             gradI = gradI &
    !             +  c_D_LEFT(ki,nk) * matmul(is_mat,rhoInkm1) + c_D_RIGHT(ki,nk) * matmul(rhoInkm1,is_mat)
    !             !   +  [sqrt(nk/abs(Ck)) * (- Ck*matmul(is_mat,rho(I_nkm1,:,:)) + conjg(Ck)*matmul(rho(I_nkm1,:,:),is_mat))]
    !             #endif
    !         endif

    !     end do
    !     ! Update the gradient array
    !     grad(I,:,:) = gradI
    ! end do
    grad(:)=0.d0

end subroutine get_gradient


end module gradient