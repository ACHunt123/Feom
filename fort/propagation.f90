! Fortran module to propagate the HEOM code
module prop_subroutines
implicit none
! Global parameters
integer(4) :: lowTcoef_switch
integer(4) :: Imax, Ktot, ns, L
real(8) :: dt, lowTcoef
real(8), parameter :: epsilon = 1.0d-12 ! epsilon for real comparisions
! Global parameter arrays
complex(8), allocatable :: C_ks(:), gam_ks(:), c_U(:,:), c_D_LEFT(:,:), c_D_RIGHT(:,:)
integer(4), allocatable :: I0s(:), lengths(:,:), ADO_index(:,:)
complex(8),allocatable :: s_mat2(:,:), is_mat(:,:), iH_mat(:,:)
! Temporary arrays for computation
complex(8), allocatable :: rhoI(:,:),rhoInkp1(:,:),rhoInkm1(:,:),gradI(:,:)
complex(8), allocatable :: k1(:,:,:), k2(:,:,:), k3(:,:,:), k4(:,:,:), ktmp(:,:,:),ADOs_tmp(:,:,:)
! Global arrays used for allocation reduction
integer(4), allocatable :: active(:), active0(:) ! active[0] is to store list of active ADOs [at start of timestep] and their indices
integer(4) :: Nactive0, Nactive ! Each vvstep, arrays are allocated (Nactive0,ns,ns), with the index of rho(I,:,:) being rho(active0(I),:,:)
logical, parameter :: prune = .false. ! if true, the code will prune ADOs with small norm
real(8), parameter :: tolerance = 1.d-20! tolerance for killing off ADOs with small norm
! Parameter for ADO printouts
logical, parameter :: print_ADOs = .false. ! if true, the code will print out the ADOs at each step
integer(4) :: nprint_ADOs = 100 ! how often to print the ADOs (every nprint_ADOs steps)
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

! Read the matrices from the files
subroutine read_matrices(ADOs)
        implicit none
        complex(8), intent(inout) :: ADOs(Imax,ns,ns)
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

subroutine ADOs_print(ADOs,Imax,ns,it)
    implicit none
    integer(4), intent(in) :: it,ns,Imax
    complex(8), intent(in) :: ADOs(Imax,ns,ns)
    real(8) :: outstr(Imax+1)
    integer(4) :: I,si
    character(len=100) :: fmt ! format string for the output
    write(fmt, '(A,I0,A)') '(E25.15,', Imax, 'E25.15)' 

    ! Open the file for writing
    open(20, file='ADOs.out', status='unknown', action='write', position='append')

    ! Write the ADOs to the file
    outstr(:) = 0.d0 ! initialize the output array
    outstr(1)=real(it)
    ! print *, 'Writing ADOs to file at timestep', it
    do I = 1, Imax
        ! do si = 1,ns
        !     outstr(I+1) = outstr(I+1) + real(ADOs(I,si,si)) ! trace calculation
        ! end do  
        outstr(I+1) = sum(abs(ADOs(I,:,:))) ! sum over all elements of the ADO
    end do
    write(20,fmt) outstr(1), outstr(2:Imax+1) ! write the ADOs to the file
    close(20)
end subroutine ADOs_print
end module prop_subroutines

program main
    use prop_subroutines
    implicit none
    complex(8), allocatable :: ADOs(:,:,:) ! we have not made this global for clarity
    ! Local variables
    real(8) :: hbar
     integer(4) :: stat,it,ki,nk,nttot,ni

    !!! LOAD PARAMETERS AND MATRICIES !!!
    ! read the parameters from the input file
    open(10, file='Fortparams', status='old', action='read', iostat=stat); read(10,*)
    read(10,'(I10, I10, D22.15, D22.15, I10, I10, D22.15, I10, I10)') Ktot, L, hbar, lowTcoef, Imax, ns, dt, nttot, lowTcoef_switch
    close(10)
    ! allocate the arrays
    allocate(iH_mat(ns,ns), is_mat(ns,ns), s_mat2(ns,ns), C_ks(Ktot), gam_ks(Ktot))
    allocate(c_U(Ktot,0:L),c_D_LEFT(Ktot,0:L),c_D_RIGHT(Ktot,0:L))
    allocate(ADO_index(Imax,Ktot), I0s(0:L+1), lengths(0:L,Ktot))
    allocate(rhoI(ns,ns),rhoInkp1(ns,ns),rhoInkm1(ns,ns),gradI(ns,ns))
    allocate(ADOs(Imax,ns,ns), k1(Imax,ns,ns), k2(Imax,ns,ns), k3(Imax,ns,ns), k4(Imax,ns,ns), ktmp(Imax,ns,ns))
    allocate(active(Imax), active0(Imax), ADOs_tmp(Imax,ns,ns))
    ! read the matrices from the files
    call read_matrices(ADOs)

    !!! FLOP REDUCTIONS !!!
    ! scale the matrices and make is2
    iH_mat = iH_mat/hbar
    is_mat = is_mat/hbar
    s_mat2 = - matmul(is_mat,is_mat) ! s_mat2 = -(i*s_mat)^2
    ! Pre-calculate the superoperator terms
    do ki = 1, Ktot
        do nk = 0, L
            c_U(ki,nk) = sqrt((nk+1)*abs(C_ks(ki)))
            if (abs(C_ks(ki))<epsilon) then! if C_ks is zero, then the superoperator term is zero (avoids 0/0 divisions)
                c_D_LEFT(ki,nk) = 0.d0 
                c_D_RIGHT(ki,nk) = 0.d0
                print *, 'Warning: C_ks(',ki,') is zero, setting superoperator terms to zero'
            else 
                c_D_LEFT(ki,nk) = -sqrt(nk/abs(C_ks(ki)))*C_ks(ki)
                c_D_RIGHT(ki,nk) = sqrt(nk/abs(C_ks(ki)))*conjg(C_ks(ki))
            end if
        end do
    end do
    ! Calculate the lengths of each block of ado indices (pascals triangle)
    do nk = 0,L 
        do ki = 1,Ktot
            lengths(nk,ki) = int(gamma(real(nk+ki-1 + 1.0D0)) / gamma(real(ki-1 + 1.0D0)) / gamma(real(nk + 1.0D0)))
        end do
    end do 
    ! Print out the hashmap for the ADO printout
    if (print_ADOs) then
        if (prune) stop 'ADO printout is not implemented for pruning'
        open(20, file='ADO_index.out', status='unknown', action='write')
        do ni = 1,Imax
            write(20,'(I10, 5I10)') ni, ADO_index(ni,:)
        end do
        close(20)
    end if

    !!! MAKE ALL ADOS INITIALLY ACTIVE !!!
    ! it has been found that setting just the first to active gives inaccurate results 
    ! this makes sense; truncation/addition of ados is done each timestep (a discrete process)
    ! if the number of ADOs explodes at t=0, then the timestep resolution is too low.
    ! also setting up this way ensures that if prune is off, the code will still work (active(I) = I)
    do ni = 1,Imax
        active(ni) = ni
        Nactive = Nactive + 1
    end do

    !!! PROPAGATION !!!
    open(10, file='output', status='unknown', action='write')
    ! Propagate the system
    do it = 1, nttot
        ! if(it==3) stop
        if (mod(it,100).eq.0) print*, it,'/',nttot, Nactive,'of',Imax,'ADOs'            ! Update the screen output
        if( mod(it,nprint_ADOs).eq.0 .and. print_ADOs) call ADOs_print(ADOs,Imax,ns,it)         ! Print the ADOs to file
        call vvstep(ADOs)
        write(10,'(5E25.15)') it*dt, &
        real(ADOs(1,1,1)), real(ADOs(1,2,2)), real(ADOs(1,1,2)), aimag(ADOs(1,1,2))
        if (abs(ADOs(1,1,1)).gt.2.d0) stop 'Density matrix has diverged'
    end do
    close(10)
    deallocate(ADOs, iH_mat, is_mat, C_ks, gam_ks, ADO_index, I0s, lengths, c_U, c_D_LEFT, c_D_RIGHT, active, s_mat2)
    deallocate(rhoI,rhoInkp1,rhoInkm1,gradI,k1,k2,k3,k4,ktmp)

end program main


