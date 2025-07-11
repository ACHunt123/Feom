program main
    use input_output
    use prop_subroutines, only: RK4step, Krylov_vecs, k1, k2, k3, k4, ktmp, ADOs_tmp, temp_grad, Krylov_dim
    use gradient, only: rhoI, rhoInkp1, rhoInkm1, gradI
    use shared_data
    implicit none
    complex(8), allocatable :: ADOs(:,:,:) ! we have not made this global for clarity
    ! Local variables
    real(8) :: hbar,time
    integer(4) :: stat,it,ki,nk,nttot,ni,ntout

    !!! LOAD PARAMETERS AND MATRICIES !!!
    ! read the parameters from the input file
    open(10, file='Fortparams', status='old', action='read', iostat=stat); read(10,*)
    read(10,'(I10, I10, D22.15, D22.15, I10, I10, D22.15, I10, I10)') Ktot, L, hbar, lowTcoef, Imax, ns, dt, nttot, lowTcoef_switch
    close(10)
    Ntot = Imax * ns * ns ! total number of elements in the ADOs array
    ! allocate the arrays
    allocate(iH_mat(ns,ns), is_mat(ns,ns), s_mat2(ns,ns), gam_ks(Ktot))
    allocate(c_U(Ktot,0:L),c_D_LEFT(Ktot,0:L),c_D_RIGHT(Ktot,0:L))
    allocate(ADO_index(Imax,Ktot), I0s(0:L+1), lengths(0:L,Ktot))
    allocate(rhoI(ns,ns),rhoInkp1(ns,ns),rhoInkm1(ns,ns),gradI(ns,ns))
    allocate(ADOs(Imax,ns,ns), k1(Imax,ns,ns), k2(Imax,ns,ns), k3(Imax,ns,ns), k4(Imax,ns,ns), ktmp(Imax,ns,ns),temp_grad(Imax,ns,ns))
    allocate(active(Imax), active0(Imax), ADOs_tmp(Imax,ns,ns))
    allocate(Krylov_vecs(0:Krylov_dim,Ntot))
    ! read the matrices from the files
    call read_matrices(ADOs)

    !!! FLOP REDUCTIONS !!!
    ! scale the matrices and make is2
    iH_mat = iH_mat/hbar
    is_mat = is_mat/hbar
    s_mat2 = - matmul(is_mat,is_mat) ! s_mat2 = -(i*s_mat)^2
    ! Calculate the lengths of each block of ado indices (pascals triangle)
    do nk = 0,L 
        do ki = 1,Ktot
            lengths(nk,ki) = int(gamma(real(nk+ki-1 + 1.0D0)) / gamma(real(ki-1 + 1.0D0)) / gamma(real(nk + 1.0D0)))
        end do
    end do 
    ! Print out the hashmap for the ADO printout
    #ifdef Print_ADOs
        #ifdef Prune
         stop 'ADO printout is not implemented for pruning'
        #endif
        open(20, file='ADO_index.out', status='unknown', action='write')
        do ni = 1,Imax
            write(20,'(I10, 5I10)') ni, ADO_index(ni,:)
        end do
        close(20)
    #endif

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
    ntout = int(max(nttot / 1000,1)) ! how often to print the output
    open(10, file='output', status='unknown', action='write')
    ! Propagate the system
    do it = 0, nttot
        if (mod(it,ntout).eq.0)  then
            time = it * dt
            print*, it,'/',nttot, Nactive,'of',Imax,'ADOs'            ! Update the screen output
            write(10,'(5E25.15)') time, &                             ! Print output to file
            real(ADOs(1,1,1)), real(ADOs(1,2,2)), real(ADOs(1,1,2)), aimag(ADOs(1,1,2))
        endif
        #ifdef Print_ADOs
        if( mod(it,nprint_ADOs).eq.0) call ADOs_print(ADOs,Imax,ns,it)         ! Print the ADOs to file
        #endif
        !verlet step
        call RK4step(ADOs)
        if (abs(ADOs(1,1,1)).gt.2.d0) stop 'Density matrix has diverged'
    end do
    close(10)
    deallocate(ADOs, iH_mat, is_mat, gam_ks, ADO_index, I0s, lengths, c_U, c_D_LEFT, c_D_RIGHT, active, s_mat2)
    deallocate(rhoI,rhoInkp1,rhoInkm1,gradI,k1,k2,k3,k4,ktmp,ADOs_tmp,active0)

end program main


