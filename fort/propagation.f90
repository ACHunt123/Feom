! Fortran module to propagate the HEOM code
module prop_subroutines
implicit none
contains


! rk4 propagation of HEOM for one step
subroutine vvstep(ADOs,ADO_index,I0s,gam_ks,C_ks,Imax,iH_mat,is_mat,Ktot,L,dt,lowTcoef,ns)
    implicit none
    integer, intent(in) :: Imax, Ktot, ns, L
    real(8), intent(in) :: dt, lowTcoef
    complex(8), intent(in) :: C_ks(Ktot), gam_ks(Ktot)
    complex(8), intent(in) :: iH_mat(ns,ns), is_mat(ns,ns)
    complex(8), intent(inout) :: ADOs(Imax,ns,ns)
    ! Arrays for the ADO index and the I0s
    integer(4), intent(in) :: ADO_index(Imax,Ktot) 
    integer(4), intent(in) :: I0s(0:L+1)
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
        complex(8) :: is_mat2(ns,ns),Ck,nk
        complex(8) :: rhoI(ns,ns),gradI(ns,ns) !temporary variables for the ADO and the gradient (speeds up code for large I and ns)
        integer(4) :: I, n_ks(Ktot), I_nkp1, I_nkm1,ki

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
                Ck = C_ks(ki)
                nk = n_ks(ki)
                ! Calculate the indices of the ADOs with nk+1 and nk-1
                I_nkp1 = I_nk_plusminus(I,ki,+1)
                I_nkm1 = I_nk_plusminus(I,ki,-1)

                if (I_nkp1.ne.-1) then
                    gradI = gradI &
                  +  sqrt((nk+1)*abs(Ck)) * ( - matmul(is_mat,rho(I_nkp1,:,:)) + matmul(rho(I_nkp1,:,:),is_mat))
                end if

                if (I_nkm1.ne.-1) then
                    gradI = gradI &
                  +  sqrt(nk/abs(Ck)) * (- Ck*matmul(is_mat,rho(I_nkm1,:,:)) + conjg(Ck)*matmul(rho(I_nkm1,:,:),is_mat))
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
        integer(4) :: indx(Ktot),I0,I1,tier0

        indx = ADO_index(I,:) !deepcopy NOT need to chekc if it is a copy and not a pointer
        tier0 = sum(indx)
        indx(ki) = indx(ki) + pm
        ! Check if the index is valid [new tier within range, no nk less than 0] - if not return -1
        if ((tier0+pm.gt.L).or.(indx(ki).lt.0)) then
            I_nk_plusminus = -1
            return
        end if
        ! Find the index in the ADOs
        ! get the block in which the index is - will speed up the search
        if (pm.eq.1) then
            I0 = I0s(tier0+1)
            I1 = I0s(tier0+2)
        else
            I0 = I0s(tier0-1)
            I1 = I0s(tier0)
        end if
        ! search for the new index within the block
        do I_nk_plusminus = I0,I1
            if (all(indx.eq.ADO_index(I_nk_plusminus,:))) return
        end do
        ! if we get this far, something is wrong
        stop 'Error in I_nk_plusminus'
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
    complex(8), allocatable :: ADOs(:,:,:), iH_mat(:,:), is_mat(:,:), C_ks(:), gam_ks(:)
    integer(4), allocatable :: ADO_index(:,:)
    integer(4), allocatable :: I0s(:)
    ! Local variables
    integer :: stat,it
    ! read the parameters from the input file
    open(10, file='Fortparams', status='old', action='read', iostat=stat); read(10,*)
    read(10,'(I10, I10, D22.15, D22.15, I10,  I10, D22.15, I10)') Ktot, L, hbar, lowTcoef, Imax, ns, dt, nttot
    close(10)
    ! print*, Ktot, L, hbar, lowTcoef, Imax, ns, dt, nttot
    ! allocate the arrays
    allocate(ADOs(Imax,ns,ns), iH_mat(ns,ns), is_mat(ns,ns), C_ks(Ktot), gam_ks(Ktot))
    allocate(ADO_index(Imax,Ktot), I0s(0:L+1))
    ! read the matrices from the files
    call read_matrices(ADOs,ADO_index,I0s,gam_ks,C_ks,Imax,iH_mat,is_mat,Ktot,L,ns)
    ! scale the matrices to reduce number of FLOPs
    iH_mat = iH_mat/hbar
    is_mat = is_mat/hbar
    lowTcoef = -lowTcoef*hbar**2 ! counterracts the scaling of the s matrix (-1 for i^2)
    ! outfile
    open(10, file='output', status='unknown', action='write')
    ! Propagate the system
    do it = 1, nttot
        if (mod(it,100).eq.0) print*, it,'/',nttot
        call vvstep(ADOs,ADO_index,I0s,gam_ks,C_ks,Imax,iH_mat,is_mat,Ktot,L,dt,lowTcoef,ns)
        write(10,*) it*dt, real(ADOs(1,1,1)), real(ADOs(1,2,2))
    end do
    close(10)
    deallocate(ADOs, iH_mat, is_mat, C_ks, gam_ks, ADO_index, I0s)
     
end program main


