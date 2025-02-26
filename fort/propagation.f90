! Fortran module to propagate the HEOM code
! module prop_subroutines
! implicit none
! contains


! rk4 propagation of HEOM for one step
subroutine vvstep(ADOs,ADO_index,I0s,gam_ks,C_ks,Imax,H_mat,s_mat,K,hbar,L,dt,lowTcoef,N_nonmats,ns)
    implicit none
    integer, intent(in) :: Imax, N_nonmats, ns, K, L
    real(8), intent(in) :: hbar,dt
    complex(8), intent(in) :: C_ks(N_nonmats+K), gam_ks(N_nonmats+K), lowTcoef
    complex(8), intent(in) :: H_mat(ns,ns), s_mat(ns,ns)
    complex(8), intent(inout) :: ADOs(Imax,ns,ns)
    ! Arrays for the ADO index and the I0s
    integer(4), intent(in) :: ADO_index(Imax,N_nonmats+K) 
    integer(4), intent(in) :: I0s(L)
    ! Local variables
    complex(8) :: ii = (0.d0,1.d0)
    complex(8) :: k1(Imax,ns,ns), k2(Imax,ns,ns), k3(Imax,ns,ns), k4(Imax,ns,ns)
    integer(4) :: nsi
    if(Imax.gt.2147483647) stop 'Imax is too large for the ADO index array'
    
    ! Calculate the k values for the Runge-Kutta method
    k1 = grad(ADOs)*dt/2.
    k2 = grad(ADOs+k1)*dt/2.
    k3 = grad(ADOs+k2)*dt
    k4 = grad(ADOs+k3)

    ! Update the density matrix
    ADOs = ADOs + (dt/6.d0)*k4 + (2.d0/3.d0)*k2 + (k3 + k1)/3.d0 ! doest work for fourth order method
    contains

    ! Function for calculating the gradient of the density, inherits the scope of the vvstep subroutine
    complex(8) function grad(rho)
        dimension grad(Imax,ns,ns)
        complex(8), intent(in) :: rho(Imax,ns,ns) ! fortran is column major so the last index is the fastest changing
        complex(8) :: s_mat2(ns,ns),Ck,nk
        integer(4) :: I, n_ks(N_nonmats+K), I_nkp1, I_nkm1,ki

        s_mat2(:,:) = matmul(s_mat,s_mat)
        grad(:,:,:) = (0.d0,0.d0)
        do I = 1, Imax
            ! Get the n values for the ADOs
            n_ks = ADO_index(I,:)

            ! Execute the diagonal superoperator terms
            grad(I,:,:) = -ii/hbar * (matmul(H_mat,rho(I,:,:)) - matmul(rho(I,:,:),H_mat)) &
                  - sum(n_ks * gam_ks) * rho(I,:,:) 
            ! Itziki Trucation (if present)
            if (int(abs(lowTcoef*10**6)).ne.0) then
                grad(I,:,:) = grad(I,:,:) - lowTcoef  &
                * (matmul(s_mat2,rho(I,:,:)) + matmul(rho(I,:,:),s_mat2) - 2.d0*matmul(matmul(s_mat,rho(I,:,:)),s_mat))
            end if

            ! Execute the off-diagonal superoperator terms
            do ki = 1, (N_nonmats+K)
                Ck = C_ks(ki)
                nk = n_ks(ki)
                ! Calculate the indices of the ADOs with nk+1 and nk-1
                I_nkp1 = I_nk_plusminus(I,ki,+1)
                I_nkm1 = I_nk_plusminus(I,ki,-1)

                if (I_nkp1.ne.-1) then
                    grad(I,:,:) = grad(I,:,:) &
                    - ii/hbar * sqrt((nk+1)*abs(Ck)) * (matmul(s_mat,rho(I_nkp1,:,:)) - matmul(rho(I_nkp1,:,:),s_mat))
                end if

                if (I_nkm1.ne.-1) then
                    grad(I,:,:) = grad(I,:,:) &
                    - ii/hbar * sqrt(nk/abs(Ck)) * (Ck*matmul(s_mat,rho(I_nkm1,:,:)) - conjg(Ck)*matmul(rho(I_nkm1,:,:),s_mat))
                end if
            end do

        end do
    end function grad


    ! Function for hashing the ADO index
    integer(4) function I_nk_plusminus(I,ki,pm) !pm is +1 or -1, I is the inital index in the ADOs, k is which nk we are changing
        implicit none
        integer(4), intent(in) :: I,ki,pm
        integer(4) :: indx(N_nonmats+K),I0,I1,tier0

        indx = ADO_index(I,:) !deepcopy NOT need to chekc if it is a copy and not a pointer
        tier0 = sum(indx)
        indx(ki) = indx(ki) + pm
        ! Check if the index is valid [new tier within range, no nk less than 0] - if not return -1
        if ((tier0+pm.gt.L).or.(indx(ki).lt.0)) then
            I_nk_plusminus = -1
            return
        end if
        ! Find the index in the ADOs
        ! get the block in which the index is
        if (pm.eq.1) then
            I0 = I0s(tier0+1)
            I1 = I0s(tier0+2)
        else
            I0 = I0s(tier0-1)
            I1 = I0s(tier0)
        end if
        ! find the index in the block
        ! print*, indx,'indx'
        ! print*, I0,I1
        do I_nk_plusminus = 1,Imax
            ! print*, ADO_index(I_nk_plusminus,:)
            if (all(indx.eq.ADO_index(I_nk_plusminus,:))) return
        end do
        ! if we get this far, something is wrong
        stop 'Error in I_nk_plusminus'
        end function I_nk_plusminus

end subroutine





! subroutine read_matrices(rho, Lv_dto2, Lk_dt, Uda, Ns, Ns2, Nx, Nx_padded)
    !     implicit none
    !     integer, intent(in) :: Ns, Ns2, Nx, Nx_padded
    !     complex(8), intent(inout) :: rho(Ns2,Nx,Nx)
    !     complex(8), intent(inout) :: Lv_dto2(Ns2,Ns2,Nx,Nx)
    !     complex(8), intent(inout) :: Lk_dt(Ns2,Ns2,Nx_padded,Nx_padded)
    !     complex(8), intent(inout) :: Uda(Ns,Ns,Nx)
    !     real(8) :: z_real, z_imag ! real and imaginary parts of the complex number
    !     integer :: xi, xj, si, sj

    !     open(20, file='Fortrho', status='old', action='read');read(20,*)
    !     open(30, file='FortLv_dto2', status='old', action='read');read(30,*)
    !     open(40, file='FortLk_dt', status='old', action='read');read(40,*)
    !     open(50, file='FortUda', status='old', action='read');read(50,*)
    !     do si = 1, Ns2 
    !     do sj = 1, Ns2
    !     do xi = 1, Nx_padded
    !     do xj = 1, Nx_padded
    !         if (si==1 .and. xi<=Nx .and. xj<=Nx) then
    !             read(20,'(D22.15)') z_real
    !             read(20,'(D22.15)') z_imag
    !             rho(sj,xi,xj) = dcmplx(z_real, z_imag)
    !         end if

    !         if (xi<=Nx .and. xj<=Nx)  then
    !             read(30,'(D22.15)') z_real
    !             read(30,'(D22.15)') z_imag
    !             Lv_dto2(si,sj,xi,xj) = dcmplx(z_real, z_imag)
    !         end if

    !         if (xj==1 .and. xi<=Nx .and. si<=Ns .and. sj<=Ns) then
    !             read(50,'(D22.15)') z_real
    !             read(50,'(D22.15)') z_imag
    !             Uda(si,sj,xi) = dcmplx(z_real, z_imag)
    !         end if
    !         read(40,'(D22.15)') z_real
    !         read(40,'(D22.15)') z_imag
    !         Lk_dt(si,sj,xi,xj) = dcmplx(z_real, z_imag)

    !     end do; end do; end do; end do  
    !     close(20);close(30);close(40);close(50) !close the files
    ! end subroutine

! end module prop_subroutines


! program main
!     use prop_subroutines
!     implicit none

     
! end program main


