module input_output
    use shared_data, only: Imax, ns, Ktot, L, hbar
    use shared_data, only: iH_mat, is_mat, gam_ks
    use shared_data, only: c_U, c_D_LEFT, c_D_RIGHT
    use shared_data, only: ADO_index, I0s, Xi_terminator
    implicit none
    contains

! Read the matrices from the files
subroutine read_matrices(ADOs)
        implicit none
        complex(8), intent(inout) :: ADOs(Imax,ns,ns)
        ! Local variables
        real(8) :: z_real, z_imag ! real and imaginary parts of the complex number
        integer :: Ii, Ij, si, sj

        ! Open the files, skipping first line
        ! small matrices
        open(11, file='Fortc_U', status='old', action='read');read(11,*)
        open(21, file='Fortc_D_LEFT', status='old', action='read');read(21,*)
        open(31, file='Fortc_D_RIGHT', status='old', action='read');read(31,*)
        open(40, file='Fortgam_ks', status='old', action='read');read(40,*)
        open(50, file='FortI0s', status='old', action='read');read(50,*)
        ! large matrices
        open(60, file='FortH_mat', status='old', action='read');read(60,*)
        open(70, file='Fortrho', status='old', action='read');read(70,*)
        open(80, file='Forts_mat', status='old', action='read');read(80,*)
        open(90, file='FortADO_index', status='old', action='read');read(90,*)
        open(100, file='FortTerminator', status='old', action='read');read(100,*)

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
                iH_mat(si,sj) = dcmplx(z_real, z_imag)*dcmplx(0.d0,1.d0)/hbar !multiply by i/hbar (as input file gives H)
                read(80,'(D22.15)') z_real
                read(80,'(D22.15)') z_imag
                is_mat(si,sj) = dcmplx(z_real, z_imag)*dcmplx(0.d0,1.d0)/hbar !multiply by i/hbar (as input file gives s)
            end if
        end do; end do; end do

        ! read the terminator superoperator
        # if LTCorr == 1
            do si = 1, ns*ns
            do sj = 1, ns*ns
                read(100,'(D22.15)') z_real
                read(100,'(D22.15)') z_imag
                Xi_terminator(1,si,sj) = dcmplx(z_real, z_imag)
            end do; end do
        # elif LTCorr == 2
            do Ii = 1, Imax
            do si = 1, ns*ns
            do sj = 1, ns*ns
                read(100,'(D22.15)') z_real
                read(100,'(D22.15)') z_imag
                Xi_terminator(Ii,si,sj) = dcmplx(z_real, z_imag)
            end do; end do; end do
        # else 
            Xi_terminator = (0.d0,0.d0)
        # endif

        ! read the small matrices
        do Ii = 1,Ktot
            read(40,'(D22.15)') z_real
            read(40,'(D22.15)') z_imag
            gam_ks(Ii) = dcmplx(z_real, z_imag)
        end do
        do Ii = 0,L+1
            read(50,'(I10)') I0s(Ii)
        end do

        !read the superoperator terms
        do Ii = 1,Ktot
            do Ij = 0,L
                read(11,'(D22.15)') z_real
                read(11,'(D22.15)') z_imag
                c_U(Ii,Ij) = dcmplx(z_real, z_imag)
                read(21,'(D22.15)') z_real
                read(21,'(D22.15)') z_imag
                c_D_LEFT(Ii,Ij) = dcmplx(z_real, z_imag)
                read(31,'(D22.15)') z_real
                read(31,'(D22.15)') z_imag
                c_D_RIGHT(Ii,Ij) = dcmplx(z_real, z_imag)
            end do
        end do
        close(11);close(21) !close the small files
        close(31);close(40);close(50);close(60);close(70);close(80);close(90) !close the files
        close(100)
    end subroutine
! print out the ADOs (if needed)
subroutine ADOs_print(ADOs,Imax,ns,it)
    implicit none
    integer(4), intent(in) :: it,ns,Imax
    complex(8), intent(in) :: ADOs(Imax,ns,ns)
    real(8) :: outstr(Imax+1)
    integer(4) :: I!,si
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
! write the system density matrix to file
subroutine writeout(funit, time, A)
    implicit none
    integer, intent(in) :: funit        ! file unit number (e.g. 10)
    real(8), intent(in) :: time         ! current time value
    complex(8), intent(in) :: A(:,:)    ! complex matrix (ns x ns)

    real(8), parameter :: tiny_cutoff = 1.0d-100
    real(8) :: reval, imval
    integer :: si, sj
    integer :: ns

    ns = size(A,1)

    ! Start the line with time
    write(funit,'(E30.15)', advance='no') time

    ! Loop through all elements of A and write on the same line
    do si = 1, ns
        do sj = 1, ns
            reval = real(A(si,sj))
            ! Clamp small values to zero to avoid exponent issues
            if (abs(reval) < tiny_cutoff) reval = 0.0d0
            ! Write number to same line
            write(funit,'(2G30.15)', advance='no') reval
        end do
    end do

    do si = 1, ns
        do sj = 1, ns
            imval = aimag(A(si,sj))
            ! Clamp small values to zero to avoid exponent issues
            if (abs(imval) < tiny_cutoff) imval = 0.0d0
            ! Write number to same line
            write(funit,'(2G30.15)', advance='no')  imval
        end do
    end do
    ! End the line (newline)
    write(funit,*)

    ! Flush the output immediately so it can be read while running
    call flush(funit)

    end subroutine
end module input_output
