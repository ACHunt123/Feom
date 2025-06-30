! Fortran module to propagate the HEOM code
module prop_subroutines
use shared_data
implicit none
! Temporary arrays for computation
complex(8), allocatable :: rhoI(:,:),rhoInkp1(:,:),rhoInkm1(:,:),gradI(:,:)
complex(8), allocatable :: k1(:,:,:), k2(:,:,:), k3(:,:,:), k4(:,:,:), ktmp(:,:,:),ADOs_tmp(:,:,:)
! Global arrays used for allocation reduction
integer(4), allocatable :: active(:), active0(:) ! active[0] is to store list of active ADOs [at start of timestep] and their indices
integer(4) :: Nactive0, Nactive ! Each vvstep, arrays are allocated (Nactive0,ns,ns), with the index of rho(I,:,:) being rho(active0(I),:,:)
contains

! rk4 propagation of HEOM for one step
subroutine vvstep(ADOs)
    implicit none
    complex(8), allocatable, intent(inout) :: ADOs(:,:,:)
    if(Imax.gt.2147483647) stop 'Imax is too large for the ADO index array'

    ! reallovate the temporary arrays if required
    if (Nactive0.ne.Nactive) then
        deallocate(k1,k2,k3,k4,ktmp)
        allocate(k1(Nactive,ns,ns), k2(Nactive,ns,ns), k3(Nactive,ns,ns), k4(Nactive,ns,ns), ktmp(Nactive,ns,ns))
    end if

    ! initialise nactive0 and active0
    ! those without 0 will be changed during propagations, those with 0 will not
    Nactive0 = Nactive
    active0 = active

    ! Calculate the k values for the Runge-Kutta method
    k1 = grad(ADOs)*dt/2. !need to split the computation here into two lines to avoid a bug
    ktmp = ADOs+k1
    k2 = grad(ktmp)*dt/2.
    ktmp = ADOs+k2
    k3 = grad(ktmp)*dt
    ktmp = ADOs+k3
    k4 = grad(ktmp)

    ! Update the density matrix
    ADOs = ADOs + (dt/6.d0)*k4 + (2.d0/3.d0)*k2 + (k3 + k1)/3.d0 ! doest work for fourth order method
    
    ! Kill off ADOs with small norm
    if (prune) call changeN(ADOs,active,active0,Nactive,Nactive0)

    contains

    ! Function for calculating the gradient of the density, inherits the scope of the vvstep subroutine
    complex(8) function grad(rho)
        dimension grad(Nactive0,ns,ns)
        complex(8), allocatable, intent(in) :: rho(:,:,:) ! fortran is column major so the last index is the fastest changing
        integer(4) :: I, n_ks(Ktot), I_nkp1, I_nkm1,ki,nk
 
        ! Loop over the ADOs
        do I = 1, Imax
            if (active0(I).eq.0) cycle !skip if not active

            ! temporary variable for the ADO (same done for the gradient, but needn't be initiallised as is overwritten)
            rhoI = rho(active0(I),:,:)
            ! Get the n values for the ADOs
            n_ks = ADO_index(I,:)

            ! Execute the diagonal superoperator terms
            gradI = ( - matmul(iH_mat,rhoI) + matmul(rhoI,iH_mat)) &
                  - sum(n_ks * gam_ks) * rhoI 

            ! Itziki Trucation (if present)
            if (lowTcoef_switch.eq.1) then
                gradI = gradI + lowTcoef  &
                * ( matmul(s_mat2,rhoI) + matmul(rhoI,s_mat2) + 2.d0*matmul(matmul(is_mat,rhoI),is_mat)) !note the + on last term is as (is)Rho(is) = - 2 s Rho s
            end if

            ! Execute the off-diagonal superoperator terms
            !NOTE: The [commented lines] are what the code is really doing (we put in the precalculated values to speed up the code)
            do ki = 1, (Ktot)
                nk = n_ks(ki)
                ! Calculate the indices of the ADOs with nk+1 and nk-1
                I_nkp1 = I_nk_plusminus(I,ki,+1)
                I_nkm1 = I_nk_plusminus(I,ki,-1)

                if (I_nkp1.ne.-1) then !check if the index is valid
                    if (active0(I_nkp1).eq.0) then      !skip if not active0 (initially active)
                        if (active(I_nkp1).eq.0) then   ! if not already active, add to the list of active ADOs that will be updated
                            Nactive = Nactive + 1
                            active(I_nkp1) = Nactive
                        end if
                    else    ! index is valid and rhos are active so do the operation
                    rhoInkp1 = rho(active0(I_nkp1),:,:)
                    gradI = gradI &
                    +  c_U(ki,nk) * ( - matmul(is_mat,rhoInkp1) + matmul(rhoInkp1,is_mat))
                    !   + [sqrt((nk+1)*abs(Ck)) * ( - matmul(is_mat,rho(I_nkp1,:,:)) + matmul(rho(I_nkp1,:,:),is_mat))]
                    endif
                endif

                if (I_nkm1.ne.-1) then
                    if(active0(I_nkm1).eq.0) then
                        if (active(I_nkm1).eq.0) then
                            Nactive = Nactive + 1
                            active(I_nkm1) = Nactive
                        end if
                    else
                    rhoInkm1 = rho(active0(I_nkm1),:,:)
                    gradI = gradI &
                    +  c_D_LEFT(ki,nk) * matmul(is_mat,rhoInkm1) + c_D_RIGHT(ki,nk) * matmul(rhoInkm1,is_mat)
                    !   +  [sqrt(nk/abs(Ck)) * (- Ck*matmul(is_mat,rho(I_nkm1,:,:)) + conjg(Ck)*matmul(rho(I_nkm1,:,:),is_mat))]
                    endif
                endif
            end do
            ! Update the gradient array
            grad(active0(I),:,:) = gradI
        end do
    end function grad


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


    ! Function for killing off ADOs that have a small norm
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

end subroutine

end module prop_subroutines