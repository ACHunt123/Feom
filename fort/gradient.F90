module gradient
use shared_data, only: Imax, ns, dt, lowTcoef_switch, lowTcoef, Nactive0, Nactive, active, active0
use shared_data, only: gam_ks, c_U, c_D_LEFT, c_D_RIGHT, ADO_index, I0s, lengths
use shared_data, only: s_mat2, is_mat, iH_mat, Ktot, L, Ntot

implicit none  
! Temporary arrays for computation
complex(8), allocatable :: rhoI(:,:),rhoInkp1(:,:),rhoInkm1(:,:),gradI(:,:)

contains

    ! Function for calculating the gradient of the density, inherits the scope of the vvstep subroutine
subroutine get_gradient(rho,grad)
    complex(8), intent(in) :: rho(Nactive0,ns,ns) ! fortran is column major so the last index is the fastest changing
    complex(8), intent(out) :: grad(Nactive0,ns,ns)
    integer(4) :: I, n_ks(Ktot), I_nkp1, I_nkm1,ki,nk

    ! Loop over the ADOs
    do I = 1, Imax
        #ifdef Prune
        if (active0(I).eq.0) cycle !skip if not active
        rhoI = rho(active0(I),:,:) ! temporary variable for the ADO (same done for the gradient, but needn't be initiallised as is overwritten)
        #else
        rhoI = rho(I,:,:)
        #endif

        ! Get the n values for the ADOs
        n_ks = ADO_index(I,:)

        ! Execute the diagonal superoperator terms
        #ifdef USEZGEMM
        call ZGEMM('N','N',ns,ns,ns,-(1.0d0,0.0d0),iH_mat,ns,rhoI,ns,(0.0d0,0.0d0),gradI,ns) ! gradI = -iH_mat * rhoI
        call ZGEMM('N','N',ns,ns,ns,(1.0d0,0.0d0),rhoI,ns,iH_mat,ns,(1.0d0,0.0d0),gradI,ns)  ! gradI = gradI + rhoI * iH_mat
        gradI = gradI  - sum(n_ks * gam_ks) * rhoI
        #else
        gradI = ( - matmul(iH_mat,rhoI) + matmul(rhoI,iH_mat)) &
                - sum(n_ks * gam_ks) * rhoI 
        #endif

        ! Itziki Trucation (if present)
        #ifdef LowTCorr
            gradI = gradI + lowTcoef  &
            * ( matmul(s_mat2,rhoI) + matmul(rhoI,s_mat2) + 2.d0*matmul(matmul(is_mat,rhoI),is_mat)) !note the + on last term is as (is)Rho(is) = - 2 s Rho s
        #endif

        ! Execute the off-diagonal superoperator terms
        !NOTE: The [commented lines] are what the code is really doing (we put in the precalculated values to speed up the code)
        do ki = 1, (Ktot)
            nk = n_ks(ki)
            ! Calculate the indices of the ADOs with nk+1 and nk-1
            I_nkp1 = I_nk_plusminus(I,ki,+1)
            I_nkm1 = I_nk_plusminus(I,ki,-1)

            if (I_nkp1.ne.-1) then !check if the index is valid
                #ifdef Prune
                if (active0(I_nkp1).eq.0) then      !skip if not active0 (initially active)
                    if (active(I_nkp1).eq.0) then   ! if not already active, add to the list of active ADOs that will be updated
                        Nactive = Nactive + 1
                        active(I_nkp1) = Nactive
                    end if
                else    ! index is valid and rhos are active so do the operation
                I_nkp1=active0(I_nkp1)
                #endif
                rhoInkp1 = rho(I_nkp1,:,:)
                #ifdef USEZGEMM
                call ZGEMM('N','N',ns,ns,ns,c_U(ki,nk),rhoInkp1,ns,is_mat,ns,(1.0d0,0.0d0),gradI,ns)    ! gradI = gradI + c_U(ki,nk) * rhoInkp1 * is_mat
                call ZGEMM('N','N',ns,ns,ns,-c_U(ki,nk),is_mat,ns,rhoInkp1,ns,(1.0d0,0.0d0),gradI,ns)   ! gradI = gradI - c_U(ki,nk) * is_mat * rhoInkp1
                #else
                gradI = gradI &
                +  c_U(ki,nk) * ( - matmul(is_mat,rhoInkp1) + matmul(rhoInkp1,is_mat))
                !   + [sqrt((nk+1)*abs(Ck)) * ( - matmul(is_mat,rho(I_nkp1,:,:)) + matmul(rho(I_nkp1,:,:),is_mat))]
                #endif
                endif
                #ifdef Prune
                endif
                #endif

            if (I_nkm1.ne.-1) then
                #ifdef Prune
                if(active0(I_nkm1).eq.0) then
                    if (active(I_nkm1).eq.0) then
                        Nactive = Nactive + 1
                        active(I_nkm1) = Nactive
                    end if
                else
                I_nkm1=active0(I_nkm1)
                #endif
                rhoInkm1 = rho(I_nkm1,:,:)
                #ifdef USEZGEMM
                call ZGEMM('N','N',ns,ns,ns,c_D_LEFT(ki,nk),is_mat,ns,rhoInkm1,ns,(1.0d0,0.0d0),gradI,ns)  ! gradI = gradI + c_D_LEFT(ki,nk) * is_mat * rhoInkm1
                call ZGEMM('N','N',ns,ns,ns,c_D_RIGHT(ki,nk),rhoInkm1,ns,is_mat,ns,(1.0d0,0.0d0),gradI,ns) ! gradI = gradI + c_D_RIGHT(ki,nk) * rhoInkm1 * is_mat
                #else
                gradI = gradI &
                +  c_D_LEFT(ki,nk) * matmul(is_mat,rhoInkm1) + c_D_RIGHT(ki,nk) * matmul(rhoInkm1,is_mat)
                !   +  [sqrt(nk/abs(Ck)) * (- Ck*matmul(is_mat,rho(I_nkm1,:,:)) + conjg(Ck)*matmul(rho(I_nkm1,:,:),is_mat))]
                #endif
                endif
                #ifdef Prune
                endif
                #endif
        end do
        ! Update the gradient array
        grad(active0(I),:,:) = gradI
    end do
    ! Recalculate the total number of elements, considering the pruning
    Ntot = Nactive*ns*ns                ! total number of elements in the ADOs array

end subroutine get_gradient


! Function for hashing the ADO index
integer(4) function I_nk_plusminus(I,ki,pm) !pm is +1 or -1, I is the inital index in the ADOs, k is which nk we are changing
    implicit none
    integer(4), intent(in) :: I,ki,pm
    integer(4) :: indx(Ktot),tier0,p,sn,np,tier1,n

    ! Get values of [n_1,n_2,...,n_Ktot] for the target ADO
    indx = ADO_index(I,:) !deepcopy NOT need to chekc if it is a copy and not a pointer
    tier0 = sum(indx)
    tier1 = tier0 + pm
    indx(ki) = indx(ki) + pm

    ! Check if the index is valid [new tier within range, no nk less than 0] - if not return -1
    if ((tier1.gt.L).or.(indx(ki).lt.0)) then
        I_nk_plusminus = -1
        return
    end if

    ! Algorithm for finding I without a hashmap          By A. C. Hunt
    !
    ! This algorithm if based on the fact that we have stored ADOs in blocks of the same leading digit, with increasing order
    ! e.g. for Ktot=3, L=2, the ADOs are stored as such:
    !
    ! [0,0,0] tier 0
    !
    ! [1,0,0]                   The First step uses I0s to find the starting index of the block (denoted by the tier)
    ! [0,1,0] tier 1
    ! [0,0,1]                   We then loop over the digits of indx (the index we are trying to find), starting left to right. 
    !                           This loop is over p =1 to Ktot-1. (no p=Ktot needed as the last digit is fixed by the tier)
    ! [2,0,0]                   
    ! [1,1,0]                   For the pth iteration, we calculate the required sum of the remaining digits required to reach the tier
    ! [1,0,1] tier 2            using the lengths array (pascals triangle) ie
    ! [0,2,0]                    for [o,x,.,.,.] {where o=correct, x=to make correct in this loop, . = remaining digits}
    ! [0,1,1]                    we calculate x = (tier - sum_k=1^p n_k), which gives the sum of the remaining digits. 
    ! [0,0,2]                    we then add the number of combinations of the remaining digits that is less than x (see the loop over n)
    !                            to I_nk_plusminus, such that the index is in the correct block --> [o,o,.,.,.]
    !                           
    !                           The output after the loop above [o,o,.,.,.] is the index of the first ADO in the block that has the leading digits [o,o]
    !                           We then repead the process for the next digit, [o,o,x,.,.] and so on until we have the correct index

    ! The algorithm is O(Ktot) in the worst case, but is usually much faster
    
    I_nk_plusminus = I0s(tier1) ! find starting index of the block
    sn = 0                      !  initialize running total of the digits 

    do p = 1,Ktot-1 ! Loop over the digit to focus on [x,.,.,.,.,.] then [o,x,.,.,.,.] etc. 
        np = indx(p) ! The value of the digit that we are looking for: n_p
        sn = sn + np ! Add this to the running total: sum_k=1^p n_k
        ! Move to the first instance where the leading digit is x = n_p
        do n = 0, tier1-sn-1    ! add number of ADOS for the sum of the remaining digits = 0,1 ... tier1-sn-1
            I_nk_plusminus = I_nk_plusminus + lengths(n,Ktot-p) ! lengths is a pascals triangle array (see main)
        end do
        if (sn==tier1) return ! If the running total of all the digits is equal to the tier, then we have found the correct index
    end do
    return
end function I_nk_plusminus

end module gradient