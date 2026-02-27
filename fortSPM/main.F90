program main
    use input_output
    use prop_subroutines, only: RK4step, k1, k2, k3, k4, ktmp!, ADOs_tmp, temp_grad
    use prop_subroutines, only: Krylov_vecs, SIAstep, Krylov_dim,L_mat,default_Krylov_dim, Recalculate_ADOs,ADOs_Krylov
    use shared_data
    use utils, only: norm
    use complex_sparse_linalg
    implicit none
    
    complex(8), allocatable :: ADOs(:) ! we have not made this global for clarity
    ! Local variables
    real(8) :: hbar,time, ADOnorm
    integer(4) :: stat,it,ki,nk,nttot,ni,ntout

    logical :: printData

    !!! LOAD PARAMETERS AND MATRICIES !!!
    ! read the parameters from the input file
    open(10, file='Fortparams.inp', status='old', action='read', iostat=stat); read(10,*)
    read(10,'(I10, D22.15, I10)') ns, dt, nttot
    close(10)
    call read_matrix('FortLiouvillian.inp', Liouvillian)
    call read_Zvec("Fortrho.inp", ADOs)
    Ntot = size(ADOs)

    
    !!! Allocations
    allocate(k1(Ntot), k2(Ntot), k3(Ntot), k4(Ntot), ktmp(Ntot))!,temp_grad(Ntot)
    
    ! Calculate the safe Krylov dimension
    Krylov_dim = min(default_Krylov_dim, Ntot)
    allocate(Krylov_vecs(Krylov_dim,Ntot), ADOs_Krylov(Krylov_dim), L_mat(Krylov_dim,Krylov_dim))           ! ADOs in the Krylov basis

    !!! Initialise SIA (such that the basis is recaluclated at the first timestep)
    ADOnorm = norm(ADOs,Ntot)
    ADOs_Krylov(:) = (1.0d0,0.0d0)*ADOnorm  ! initialise the ADOs in the Krylov basis (made to trigger recalculation initially)
    Krylov_vecs(:,:) = (0.0d0,0.0d0)        ! initialise the Krylov vectors
    Krylov_vecs(1,:) = ADOs/ADOnorm ! set the first Krylov vector top be the ADOs


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
            call writeout(10,time,ADOs(1:ns**2))
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
        ! if (abs(ADOs(1)).gt.2.d0) stop 'Density matrix has diverged'
    end do

    close(10)
    deallocate(ADOs)
    ! need to deallocate everything 
    ! deallocate(k1,k2,k3,k4,ktmp)

end program main


