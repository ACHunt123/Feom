module input_output
    implicit none
    contains

! Read a vectorized complex matrix from a file
subroutine read_Zvec(filename, vec)
        implicit none
        ! Arguments
        character(len=*), intent(in) :: filename  
        complex(8), allocatable, intent(out)      :: vec(:)
        ! Local variables
        integer :: n_total   
        real(8) :: z_real, z_imag 
        integer :: i
        ! Open the file
        open(70, file=filename, status='old', action='read')
        read(70, *) ! Skip the first line (header)
        read(70, *) n_total ! Read the length
        allocate(vec(n_total))
        ! Read the data
        do i = 1, n_total
            read(70, '(D22.15)') z_real
            read(70, '(D22.15)') z_imag
            vec(i) = dcmplx(z_real, z_imag)
        end do
        ! Close the file
        close(70)
    end subroutine read_Zvec

subroutine ADOs_print(ADOs, ns, Ntot, time)
    implicit none
    integer(4), intent(in) :: ns,Ntot
    real(8), intent(in) :: time
    complex(8), intent(in) :: ADOs(Ntot) 
    integer(4) :: n, Nados, start_idx, end_idx
    integer(4), parameter :: funit = 11
    
    Nados = Ntot / (ns**2) ! get the number of ADOs
    open(unit=funit, file='ADOs.dat', status='unknown', action='write', position='append')
    ! print out each ADO sequentially
    do n = 1, Nados
        start_idx = (n-1) * (ns**2) + 1
        end_idx   = n * (ns**2)
        call writeout(funit, time, ADOs(start_idx:end_idx))
    end do
    close(funit)
    end subroutine ADOs_print

! write the system density matrix to file
subroutine writeout(funit, time, A)
    implicit none
    integer, intent(in) :: funit        ! file unit number (e.g. 10)
    real(8), intent(in) :: time         ! current time value
    complex(8), intent(in) :: A(:)      ! complex 1D array of unknown length

    real(8), parameter :: tiny_cutoff = 1.0d-100
    real(8) :: reval, imval
    integer :: i
    integer :: n

    ! Get the length of the 1D array
    n = size(A)
    ! Start the line with time
    write(funit,'(E30.15E3)', advance='no') time
    ! Loop through all elements of A and write real parts on the same line
    do i = 1, n
        reval = real(A(i))
        ! Clamp small values to zero to avoid exponent issues
        if (abs(reval) < tiny_cutoff) reval = 0.0d0
        ! Write number to same line (Fixed format from 2G to G)
        write(funit,'(G30.15E3)', advance='no') reval
    end do
    ! Loop through all elements of A and write imaginary parts on the same line
    do i = 1, n
        imval = aimag(A(i))
        ! Clamp small values to zero to avoid exponent issues
        if (abs(imval) < tiny_cutoff) imval = 0.0d0
        ! Write number to same line (Fixed format from 2G to G)
        write(funit,'(G30.15E3)', advance='no') imval
    end do
    ! End the line (newline)
    write(funit,*)
    ! Flush the output immediately so it can be read while running
    call flush(funit)
end subroutine writeout


end module input_output
