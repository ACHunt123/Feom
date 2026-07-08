program main
    use input_output
    use integrator, only: step_forward,init_integrator,cleanup_integrator,Recalculate_ADOs
    use shared_data
    use utils, only: norm
    use complex_sparse_linalg, only: destroy_matrix
    implicit none
    
    complex(8), allocatable :: ADOs(:) ! we have not made this global for clarity
    ! Local variables
    real(8) :: time, ADOnorm
    integer(4) :: stat,it,nttot,ntout
    logical :: printData

    !!! LOAD PARAMETERS AND MATRICIES !!!
    open(10, file='Fortparams.inp', status='old', action='read', iostat=stat); read(10,*)
    read(10,'(I10, D22.15, I10, I10)') ns, dt, nttot, ntout
    close(10)
    call read_matrix('FortLiouvillian.inp', Liouvillian)
    call read_Zvec("Fortrho.inp", ADOs)
    Ntot = size(ADOs)

    !!! Allocations for work matrices
    call init_integrator(ADOs)

    !!! PROPAGATION !!!
    nttot = (nttot / ntout) * ntout  ! force the last timestep to be printed out
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
        time = it * dt
        if( mod(it,nprint_ADOs).eq.0) call ADOs_print(ADOs,ns,Ntot,time)         ! Print the ADOs to file
#endif

        !!! verlet step (skipping last one)
        if (it < nttot) then
            call step_forward(ADOs)
        endif
        !!! Check for divergence
        if (abs(ADOs(1)).gt.2.d0) stop 'Density matrix has diverged'
    end do

    close(10)
    !!! Write final state to file
    call write_Zvec("Fortrho.oup", ADOs)

    !!! Deallocations 
    deallocate(ADOs)
    call destroy_matrix(Liouvillian)

    call cleanup_integrator()
end program main


