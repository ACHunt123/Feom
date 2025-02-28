! Fortran module to propagate the HEOM code
module prop_subroutines
implicit none
contains


! rk4 propagation of HEOM for one step
subroutine vvstep(ADOs,ADO_index,I0s,lengths,gam_ks,c_U,C_D_LEFT,c_D_RIGHT,Imax,iH_mat,is_mat,Ktot,L,dt,lowTcoef,ns)
    implicit none
    integer, intent(in) :: Imax, Ktot, ns, L
    real(8), intent(in) :: dt, lowTcoef
    complex(8), intent(in) :: gam_ks(Ktot), c_U(Ktot,0:L),c_D_LEFT(Ktot,0:L),c_D_RIGHT(Ktot,0:L)   
    complex(8), intent(in) :: iH_mat(ns,ns), is_mat(ns,ns)
    complex(8), intent(inout) :: ADOs(Imax,ns,ns)
    ! Arrays for the ADO index and the I0s
    integer(4), intent(in) :: ADO_index(Imax,Ktot) 
    integer(4), intent(in) :: I0s(0:L+1), lengths(0:L,Ktot)
    ! Local variables
    complex(8) :: k1(Imax,ns,ns), k2(Imax,ns,ns), k3(Imax,ns,ns), k4(Imax,ns,ns), ktmp(Imax,ns,ns)
    if(Imax.gt.2147483647) stop 'Imax is too large for the ADO index array'
    
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
    contains

    ! Function for calculating the gradient of the density, inherits the scope of the vvstep subroutine
    complex(8) function grad(rho)
        dimension grad(Imax,ns,ns)
        complex(8), intent(in) :: rho(Imax,ns,ns) ! fortran is column major so the last index is the fastest changing
        complex(8) :: is_mat2(ns,ns)!,Ck!,nk
        complex(8) :: rhoI(ns,ns),gradI(ns,ns), rhoInkp1(ns,ns),rhoInkm1(ns,ns) !temporary variables for the ADO and the gradient (speeds up code for large I and ns)
        integer(4) :: I, n_ks(Ktot), I_nkp1, I_nkm1,ki,nk
 
        is_mat2(:,:) = matmul(is_mat,is_mat)

        ! Loop over the ADOs
        do I = 1, Imax
            ! temporary variable for the ADO (same done for the gradient, but needn't be initiallised as is overwritten)
            rhoI = rho(I,:,:)
            ! Get the n values for the ADOs
            n_ks = ADO_index(I,:)

            ! Execute the diagonal superoperator terms
            gradI = ( - matmul(iH_mat,rhoI) + matmul(rhoI,iH_mat)) &
                  - sum(n_ks * gam_ks) * rhoI 

            ! Itziki Trucation (if present)
            if (int(abs(lowTcoef*10**6)).ne.0) then
                gradI = gradI + lowTcoef  &
                * (- matmul(is_mat2,rhoI) - matmul(rhoI,is_mat2) + 2.d0*matmul(matmul(is_mat,rhoI),is_mat))
            end if

            ! Execute the off-diagonal superoperator terms
            do ki = 1, (Ktot)
                nk = n_ks(ki)
                ! Calculate the indices of the ADOs with nk+1 and nk-1
                I_nkp1 = I_nk_plusminus(I,ki,+1)
                I_nkm1 = I_nk_plusminus(I,ki,-1)

                !NOTE: The commented lines are what the code is really doing (we put in the precalculated values to speed up the code)
                if (I_nkp1.ne.-1) then
                    rhoInkp1 = rho(I_nkp1,:,:)
                    gradI = gradI &
                    +  c_U(ki,nk) * ( - matmul(is_mat,rhoInkp1) + matmul(rhoInkp1,is_mat))
                    !   +  sqrt((nk+1)*abs(Ck)) * ( - matmul(is_mat,rho(I_nkp1,:,:)) + matmul(rho(I_nkp1,:,:),is_mat))
                end if

                if (I_nkm1.ne.-1) then
                    rhoInkm1 = rho(I_nkm1,:,:)
                    gradI = gradI &
                !   +  sqrt(nk/abs(Ck)) * (- Ck*matmul(is_mat,rho(I_nkm1,:,:)) + conjg(Ck)*matmul(rho(I_nkm1,:,:),is_mat))
                  +  c_D_LEFT(ki,nk) * matmul(is_mat,rhoInkm1) + c_D_RIGHT(ki,nk) * matmul(rhoInkm1,is_mat)
                end if
            end do
            ! Update the gradient array
            grad(I,:,:) = gradI
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

end subroutine

! Read the matrices from the files
subroutine read_matrices(ADOs,ADO_index,I0s,gam_ks,C_ks,Imax,iH_mat,is_mat,Ktot,L,ns)
        implicit none
        integer, intent(in) :: Imax, ns, Ktot, L
        complex(8), intent(inout) :: C_ks(Ktot), gam_ks(Ktot)
        complex(8), intent(inout) :: iH_mat(ns,ns), is_mat(ns,ns)
        complex(8), intent(inout) :: ADOs(Imax,ns,ns)
        ! Arrays for the ADO index and the I0s
        integer(4), intent(inout) :: ADO_index(Imax,Ktot) 
        integer(4), intent(inout) :: I0s(0:L+1)
        ! Local variables
        real(8) :: z_real, z_imag ! real and imaginary parts of the complex number
        integer :: Ii, Ij, si, sj

        ! Open the files, skipping first line
        ! small matrices
        open(30, file='FortC_ks', status='old', action='read');read(30,*)
        open(40, file='Fortgam_ks', status='old', action='read');read(40,*)
        open(50, file='FortI0s', status='old', action='read');read(50,*)
        ! large matrices
        open(60, file='FortH_mat', status='old', action='read');read(60,*)
        open(70, file='Fortrho', status='old', action='read');read(70,*)
        open(80, file='Forts_mat', status='old', action='read');read(80,*)
        open(90, file='FortADO_index', status='old', action='read');read(90,*)

        ! read the large matrices
        do Ii = 1, Imax
            do Ij = 1, Ktot
                read(90,'(I10)') ADO_index(Ii,Ij)
            end do
        do si = 1, ns
        do sj = 1, ns
            read(70,'(D22.15)') z_real
            read(70,'(D22.15)') z_imag
            ADOs(Ii,si,sj) = dcmplx(z_real, z_imag)
            if (Ii==1) then
                read(60,'(D22.15)') z_real
                read(60,'(D22.15)') z_imag
                iH_mat(si,sj) = dcmplx(z_real, z_imag)*dcmplx(0.d0,1.d0) !multiply by i (as input file gives H)
                read(80,'(D22.15)') z_real
                read(80,'(D22.15)') z_imag
                is_mat(si,sj) = dcmplx(z_real, z_imag)*dcmplx(0.d0,1.d0) !multiply by i (as input file gives s)
            end if
        end do; end do; end do

        ! read the small matrices
        do Ii = 1,Ktot
            read(30,'(D22.15)') z_real
            read(30,'(D22.15)') z_imag
            C_ks(Ii) = dcmplx(z_real, z_imag)
            read(40,'(D22.15)') z_real
            read(40,'(D22.15)') z_imag
            gam_ks(Ii) = dcmplx(z_real, z_imag)
        end do
        do Ii = 0,L+1
            read(50,'(I10)') I0s(Ii)
        end do

        close(30);close(40);close(50);close(60);close(70);close(80);close(90) !close the files
    end subroutine

end module prop_subroutines

program main
    use prop_subroutines
    implicit none
    integer :: Imax, Ktot, ns, L, nttot
    real(8) :: hbar, dt, lowTcoef
    complex(8), allocatable :: ADOs(:,:,:), iH_mat(:,:), is_mat(:,:), C_ks(:), gam_ks(:), c_U(:,:), c_D_LEFT(:,:), c_D_RIGHT(:,:)
    integer(4), allocatable :: ADO_index(:,:)
    integer(4), allocatable :: I0s(:), lengths(:,:)
    ! Local variables
    integer :: stat,it,ki,nk
    ! read the parameters from the input file
    open(10, file='Fortparams', status='old', action='read', iostat=stat); read(10,*)
    read(10,'(I10, I10, D22.15, D22.15, I10,  I10, D22.15, I10)') Ktot, L, hbar, lowTcoef, Imax, ns, dt, nttot
    close(10)
    ! print*, Ktot, L, hbar, lowTcoef, Imax, ns, dt, nttot
    ! allocate the arrays
    allocate(ADOs(Imax,ns,ns), iH_mat(ns,ns), is_mat(ns,ns), C_ks(Ktot), gam_ks(Ktot))
    allocate(c_U(Ktot,0:L),c_D_LEFT(Ktot,0:L),c_D_RIGHT(Ktot,0:L))
    allocate(ADO_index(Imax,Ktot), I0s(0:L+1), lengths(0:L,Ktot))
    ! read the matrices from the files
    call read_matrices(ADOs,ADO_index,I0s,gam_ks,C_ks,Imax,iH_mat,is_mat,Ktot,L,ns)
    ! FLOP REDUCTIONS
    ! scale the matrices
    iH_mat = iH_mat/hbar
    is_mat = is_mat/hbar
    lowTcoef = -lowTcoef*hbar**2 ! counterracts the scaling of the s matrix (-1 for i^2)
    ! Pre-calculate the superoperator terms
    do ki = 1, Ktot
        do nk = 0, L
            c_U(ki,nk) = sqrt((nk+1)*abs(C_ks(ki)))
            c_D_LEFT(ki,nk) = -sqrt(nk/abs(C_ks(ki)))*C_ks(ki)
            c_D_RIGHT(ki,nk) = sqrt(nk/abs(C_ks(ki)))*conjg(C_ks(ki))
        end do
    end do
    ! Calculate the lengths of each block of ado indices (pascals triangle)
    do nk = 0,L 
        do ki = 1,Ktot
            lengths(nk,ki) = int(gamma(real(nk+ki-1 + 1.0D0)) / gamma(real(ki-1 + 1.0D0)) / gamma(real(nk + 1.0D0)))
        end do
    end do 

    ! Ready to go
    open(10, file='output', status='unknown', action='write')
    ! Propagate the system
    do it = 1, nttot
        if (mod(it,100).eq.0) print*, it,'/',nttot
        call vvstep(ADOs,ADO_index,I0s,lengths,gam_ks,c_U,C_D_LEFT,c_D_RIGHT,Imax,iH_mat,is_mat,Ktot,L,dt,lowTcoef,ns)
        write(10,'(5E25.15)') it*dt, &
        real(ADOs(1,1,1)), real(ADOs(1,2,2)), real(ADOs(1,1,2)), aimag(ADOs(1,1,2))
        if (abs(ADOs(1,1,1)).gt.1.d8) stop 'Density matrix has diverged'
    end do
    close(10)
    deallocate(ADOs, iH_mat, is_mat, C_ks, gam_ks, ADO_index, I0s, lengths, c_U, c_D_LEFT, c_D_RIGHT)
     
end program main


