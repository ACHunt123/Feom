program main
    use input_output
    use prop_subroutines
    use shared_data
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


