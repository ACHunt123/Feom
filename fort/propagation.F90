program main
    use input_output
    use prop_subroutines, only: RK4step, k1, k2, k3, k4, ktmp, ADOs_tmp, temp_grad
    use prop_subroutines, only: Krylov_vecs, SIAstep, Krylov_dim, Recalculate_ADOs,ADOs_Krylov
    use gradient, only: rhoI, rhoInkp1, rhoInkm1, gradI, result_vec, rhoI_vec
    use shared_data
    use utils, only: norm
    implicit none
    complex(8), allocatable :: ADOs(:,:,:) ! we have not made this global for clarity
    ! Local variables
    real(8) :: time, ADOnorm
    integer(4) :: stat,it,ki,nk,nttot,ni,ntout
    logical :: printData

    !!! LOAD PARAMETERS AND MATRICIES !!!
    ! read the parameters from the input file
    open(10, file='Fortparams', status='old', action='read', iostat=stat); read(10,*)
    ! read(10,'(I10, I10, D22.15, I10, I10, D22.15, I10)') Ktot, L, hbar, Imax, ns, dt, nttot ! old FIXED WIDTH FORMAT
    read(10,*) Ktot, L, hbar, Imax, ns, dt, nttot
    close(10)
    Ntot = Imax * ns * ns ! total number of elements in the ADOs array
    ! allocate the arrays
    allocate(iH_mat(ns,ns), is_mat(ns,ns), gam_ks(Ktot))
    allocate(c_U(Ktot,0:L),c_D_LEFT(Ktot,0:L),c_D_RIGHT(Ktot,0:L))
    allocate(ADO_index(Imax,Ktot), I0s(0:L+1), lengths(0:L,Ktot))
    allocate(rhoI(ns,ns),rhoInkp1(ns,ns),rhoInkm1(ns,ns),gradI(ns,ns),result_vec(ns*ns),rhoI_vec(ns*ns))
    allocate(ADOs(Imax,ns,ns), k1(Imax,ns,ns), k2(Imax,ns,ns), k3(Imax,ns,ns), k4(Imax,ns,ns), ktmp(Imax,ns,ns),temp_grad(Imax,ns,ns))
    allocate(ADOs_tmp(Imax,ns,ns))
    allocate(Krylov_vecs(Krylov_dim,Ntot))
    # if LTCorr == 0
        allocate(Xi_terminator(1,1,1))
    # elif LTCorr == 1
        allocate(Xi_terminator(1,ns*ns,ns*ns))
    # elif LTCorr == 2
        allocate(Xi_terminator(Imax,ns*ns,ns*ns))
    # endif
    ! read the matrices from the files
    call read_matrices(ADOs)

    !!! FLOP REDUCTIONS AND PRINTOUTS !!!

    ! Calculate the lengths of each block of ado indices (pascals triangle)
    do nk = 0,L 
        do ki = 1,Ktot
            lengths(nk,ki) = int(gamma(real(nk+ki-1 + 1.0D0)) / gamma(real(ki-1 + 1.0D0)) / gamma(real(nk + 1.0D0)))
        end do
    end do 
    ! Print out the hashmap for the ADO printout
    #ifdef Print_ADOs
        open(20, file='ADO_index.out', status='unknown', action='write')
        do ni = 1,Imax
            write(20,'(I10, 5I10)') ni, ADO_index(ni,:)
        end do
        close(20)
    #endif


    !!! Initialise SIA (such that the basis is recaluclated at the first timestep)
    ADOnorm = norm(ADOs,Ntot)
    ADOs_Krylov(:) = (1.0d0,0.0d0)*ADOnorm  ! initialise the ADOs in the Krylov basis (made to trigger recalculation initially)
    Krylov_vecs(:,:) = (0.0d0,0.0d0)        ! initialise the Krylov vectors
    Krylov_vecs(1,:) = reshape(ADOs(:,:,:), [Ntot])/ADOnorm ! set the first Krylov vector top be the ADOs

    !!! PROPAGATION !!!
    ntout = int(max(nttot / 1000,1)) ! how often to print the output
    open(10, file='output', status='unknown', action='write')
    ! Propagate the system
    do it = 0, nttot
        printData = mod(it,ntout).eq.0
        if (printData)  then
            time = it * dt
            #ifdef SIA 
            call Recalculate_ADOs(ADOs) ! Recalculate the ADOs from the Krylov subspace
            #endif
            
            !!! Write output and message to screen
            print*, it,'/',nttot           ! Update the screen output
            ! write(10,'(5E25.15)') time, &                             ! Print output to file
            ! real(ADOs(1,1,1)), real(ADOs(1,2,2)), real(ADOs(1,1,2)), aimag(ADOs(1,1,2))
            write(10,'(5000E25.15)') time, &
                real(reshape(ADOS(1,:,:), [ns**2])), &
                aimag(reshape(ADOS(1,:,:), [ns**2]))
        endif
        #ifdef Print_ADOs
        if( mod(it,nprint_ADOs).eq.0) call ADOs_print(ADOs,Imax,ns,it)         ! Print the ADOs to file
        #endif

        !!! verlet step
        #ifdef SIA
        call SIAstep(ADOs)
        #else
        call RK4step(ADOs)
        #endif

        !!! Check for divergence
        if (abs(ADOs(1,1,1)).gt.2.d0) stop 'Density matrix has diverged'
    end do
    close(10)
    
    deallocate(ADOs, iH_mat, is_mat, gam_ks, ADO_index, I0s, lengths, c_U, c_D_LEFT, c_D_RIGHT)
    deallocate(rhoI,rhoInkp1,rhoInkm1,gradI,k1,k2,k3,k4,ktmp,ADOs_tmp)
    deallocate(temp_grad,Krylov_vecs,Xi_terminator)
    call exit(0)
    
end program main


