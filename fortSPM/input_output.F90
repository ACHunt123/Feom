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
end module input_output
