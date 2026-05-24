module precision_module
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 300)
end module precision_module

module constants_module
  use precision_module
  implicit none

  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp

  ! hbar^2/(2m_N) in MeV fm^2.
  real(dp), parameter :: hbar2_over_2m = 20.73553_dp

  ! Simplified SkP parameters used in the Supplement of Wu et al.
  ! Units: t0 [MeV fm^3], t1,t2 [MeV fm^5], t3 [MeV fm^(3+3 gamma)].
  real(dp), parameter :: t0_skp = -2931.70_dp
  real(dp), parameter :: t1_skp =   320.62_dp
  real(dp), parameter :: t2_skp =  -337.41_dp
  real(dp), parameter :: t3_skp = 18708.97_dp
  real(dp), parameter :: x2_skp =    -0.53732_dp
  real(dp), parameter :: gamma_skp = 1.0_dp/6.0_dp

  ! KEDF constants used in the non-local NLF postprocessing.
  ! tau_TF = c_TF rho^(5/3), where T_TF = hbar^2/(2m) int tau_TF d^3r.
  real(dp), parameter :: c_TF = (3.0_dp/5.0_dp) * (3.0_dp*pi*pi)**(2.0_dp/3.0_dp)

  ! Eq. (21)-type coefficients, without hbar^2/(2m). The suffix NL denotes
  ! the nonlocal two-point term; CORR denotes the subtractive local correction
  ! that restores the Thomas-Fermi uniform-matter limit.
  ! Kernel-strength rescaling (dd-consistent recalibration):
  !   tau_lambda = tau_vW + lambda c_NL0 rho^(5/6) (f1*rho^(5/6))
  !                       - (lambda c_NL0 - c_TF) rho^(5/3),
  ! which preserves the TF (uniform-matter) limit tau -> c_TF rho^(5/3)
  ! exactly for every lambda; lambda = 1 restores the original coefficients.
  real(dp), parameter :: c_NL0   = (96.0_dp/125.0_dp) * (3.0_dp*pi*pi)**(2.0_dp/3.0_dp)
  real(dp), parameter :: c_CORR0 = (21.0_dp/125.0_dp) * (3.0_dp*pi*pi)**(2.0_dp/3.0_dp)
  real(dp) :: c_NL   = c_NL0
  real(dp) :: c_CORR = c_CORR0
  real(dp) :: kernel_lambda_saved = 1.0_dp
  ! Spatial-gradient regularization coefficient for the shared multi-nucleus fit:
  !   E_grad = s2_grad * (hbar^2/2m) * int rho^(1/3) |grad rho|^2 d^3r.
  ! Vanishes for uniform matter (grad rho = 0), so it does not change the
  ! uniform-matter (Thomas-Fermi) limit; the rho^(1/3) prefactor keeps the
  ! Euler-Lagrange potential finite in the low-density tail. Default 0.
  real(dp) :: s2_grad = 0.0_dp

  real(dp), parameter :: rho_floor  = 1.0e-20_dp
  real(dp), parameter :: tiny_denom = 1.0e-30_dp

contains
  subroutine set_kernel_lambda(lam)
    real(dp), intent(in) :: lam
    c_NL   = lam * c_NL0
    c_CORR = c_NL - c_TF          ! keeps c_NL - c_CORR = c_TF (TF limit)
    kernel_lambda_saved = lam
    if (abs(lam - 1.0_dp) > 1.0e-14_dp) then
       write(*,'(A,F10.5,A,ES14.6,A,ES14.6)') &
            'Kernel rescaling lambda = ', lam, '  c_NL = ', c_NL, '  c_CORR = ', c_CORR
    end if
  end subroutine set_kernel_lambda
end module constants_module

module radial_grid_module
  use precision_module
  use constants_module, only: pi
  implicit none
contains

  subroutine make_uniform_interior_grid(nr, rbox, r, dr, w)
    integer, intent(in) :: nr
    real(dp), intent(in) :: rbox
    real(dp), intent(out) :: r(nr), dr, w(nr)
    integer :: i

    if (nr < 3) stop 'ERROR: nr must be >= 3.'
    if (rbox <= 0.0_dp) stop 'ERROR: rbox must be positive.'

    ! Dirichlet boundaries are at r=0 and r=rbox.
    ! The unknowns u_i live on interior points i=1,...,nr.
    dr = rbox / real(nr + 1, dp)
    do i = 1, nr
       r(i) = real(i, dp) * dr
       w(i) = dr
    end do
  end subroutine make_uniform_interior_grid


  real(dp) function integrate_r3(n, r, w, f)
    integer, intent(in) :: n
    real(dp), intent(in) :: r(n), w(n), f(n)
    integer :: i

    integrate_r3 = 0.0_dp
    do i = 1, n
       integrate_r3 = integrate_r3 + 4.0_dp*pi*r(i)**2 * f(i) * w(i)
    end do
  end function integrate_r3


  real(dp) function integrate_1d(n, w, f)
    integer, intent(in) :: n
    real(dp), intent(in) :: w(n), f(n)
    integer :: i

    integrate_1d = 0.0_dp
    do i = 1, n
       integrate_1d = integrate_1d + f(i) * w(i)
    end do
  end function integrate_1d


  subroutine first_derivative_uniform(n, dr, f, df)
    integer, intent(in) :: n
    real(dp), intent(in) :: dr, f(n)
    real(dp), intent(out) :: df(n)
    integer :: i

    df(1) = (f(2) - f(1)) / dr
    do i = 2, n-1
       df(i) = (f(i+1) - f(i-1)) / (2.0_dp*dr)
    end do
    df(n) = (f(n) - f(n-1)) / dr
  end subroutine first_derivative_uniform


  subroutine second_derivative_uniform(n, dr, f, d2f)
    integer, intent(in) :: n
    real(dp), intent(in) :: dr, f(n)
    real(dp), intent(out) :: d2f(n)
    integer :: i

    ! Simple one-sided boundary approximations.
    d2f(1) = (f(3) - 2.0_dp*f(2) + f(1)) / (dr*dr)
    do i = 2, n-1
       d2f(i) = (f(i+1) - 2.0_dp*f(i) + f(i-1)) / (dr*dr)
    end do
    d2f(n) = (f(n) - 2.0_dp*f(n-1) + f(n-2)) / (dr*dr)
  end subroutine second_derivative_uniform


  subroutine radial_laplacian_uniform(n, r, dr, f, lap)
    integer, intent(in) :: n
    real(dp), intent(in) :: r(n), dr, f(n)
    real(dp), intent(out) :: lap(n)

    real(dp), allocatable :: df(:), d2f(:)
    integer :: i

    allocate(df(n), d2f(n))
    call first_derivative_uniform(n, dr, f, df)
    call second_derivative_uniform(n, dr, f, d2f)

    do i = 1, n
       lap(i) = d2f(i) + 2.0_dp * df(i) / max(r(i), 1.0e-10_dp)
    end do

    deallocate(df, d2f)
  end subroutine radial_laplacian_uniform


  subroutine compute_tf_vw_tau(n, r, dr, rho, tau_tf, tau_vw)
    use constants_module, only: c_TF, rho_floor
    integer, intent(in) :: n
    real(dp), intent(in) :: r(n), dr, rho(n)
    real(dp), intent(out) :: tau_tf(n), tau_vw(n)

    real(dp), allocatable :: sqrtrho(:), dsqrtrho(:)
    integer :: i

    allocate(sqrtrho(n), dsqrtrho(n))

    do i = 1, n
       sqrtrho(i) = sqrt(max(rho(i), rho_floor))
       tau_tf(i)  = c_TF * max(rho(i), 0.0_dp)**(5.0_dp/3.0_dp)
    end do

    call first_derivative_uniform(n, dr, sqrtrho, dsqrtrho)

    do i = 1, n
       tau_vw(i) = dsqrtrho(i)**2
    end do

    deallocate(sqrtrho, dsqrtrho)
  end subroutine compute_tf_vw_tau

end module radial_grid_module

module skyrme_skp_module
  use precision_module
  use constants_module
  use radial_grid_module, only: radial_laplacian_uniform
  implicit none
contains

  subroutine compute_skp_potential(nr, r, dr, rho_tot, vks)
    integer, intent(in) :: nr
    real(dp), intent(in) :: r(nr), dr, rho_tot(nr)
    real(dp), intent(out) :: vks(nr)

    real(dp), allocatable :: lap(:)
    real(dp) :: cgrad, rh
    integer :: i

    allocate(lap(nr))
    call radial_laplacian_uniform(nr, r, dr, rho_tot, lap)

    cgrad = 9.0_dp*t1_skp - 5.0_dp*t2_skp - 4.0_dp*t2_skp*x2_skp

    do i = 1, nr
       rh = max(rho_tot(i), 0.0_dp)
       vks(i) = 0.75_dp*t0_skp*rh &
            + ((2.0_dp + gamma_skp)/16.0_dp)*t3_skp*rh**(1.0_dp + gamma_skp) &
            - (1.0_dp/32.0_dp)*cgrad*lap(i)
    end do

    deallocate(lap)
  end subroutine compute_skp_potential

end module skyrme_skp_module

module lapack_tridiag_module
  use precision_module
  implicit none

  interface
     subroutine dstev(jobz, n, d, e, z, ldz, work, info)
       import :: dp
       character(len=1), intent(in) :: jobz
       integer, intent(in) :: n, ldz
       real(dp), intent(inout) :: d(*)
       real(dp), intent(inout) :: e(*)
       real(dp), intent(out) :: z(ldz,*)
       real(dp), intent(out) :: work(*)
       integer, intent(out) :: info
     end subroutine dstev
  end interface

contains

  subroutine diagonalize_tridiagonal(n, diag, offdiag, eigval, eigvec)
    integer, intent(in) :: n
    real(dp), intent(in) :: diag(n), offdiag(n-1)
    real(dp), intent(out) :: eigval(n), eigvec(n,n)

    real(dp), allocatable :: d(:), e(:), work(:)
    integer :: info

    allocate(d(n), e(n-1), work(max(1,2*n-2)))
    d(:) = diag(:)
    e(:) = offdiag(:)

    call dstev('V', n, d, e, eigvec, n, work, info)
    if (info /= 0) then
       write(*,*) 'ERROR: dstev failed with INFO = ', info
       stop
    end if

    eigval(:) = d(:)

    deallocate(d, e, work)
  end subroutine diagonalize_tridiagonal

end module lapack_tridiag_module

module ks_orbital_module
  use precision_module
  use constants_module
  use lapack_tridiag_module
  implicit none
contains

  subroutine solve_radial_l(nr, r, dr, l, vks, eigval, uvec)
    integer, intent(in) :: nr, l
    real(dp), intent(in) :: r(nr), dr, vks(nr)
    real(dp), intent(out) :: eigval(nr), uvec(nr,nr)

    real(dp), allocatable :: diag(:), offdiag(:), z(:,:)
    real(dp) :: c, rr
    integer :: i

    allocate(diag(nr), offdiag(nr-1), z(nr,nr))

    c = hbar2_over_2m

    do i = 1, nr
       rr = max(r(i), 1.0e-10_dp)
       diag(i) = 2.0_dp*c/(dr*dr) &
            + c*real(l*(l+1), dp)/(rr*rr) &
            + vks(i)
    end do

    do i = 1, nr-1
       offdiag(i) = -c/(dr*dr)
    end do

    call diagonalize_tridiagonal(nr, diag, offdiag, eigval, z)

    ! dstev returns eigenvectors normalized as sum_i z_i^2 = 1.
    ! We use radial wave functions normalized as int u(r)^2 dr = 1.
    uvec(:,:) = z(:,:) / sqrt(dr)

    deallocate(diag, offdiag, z)
  end subroutine solve_radial_l

end module ks_orbital_module

module sorting_module
  use precision_module
  implicit none
contains

  subroutine sort_levels(nlev, energy, lval, idx_l_state)
    integer, intent(in) :: nlev
    real(dp), intent(inout) :: energy(nlev)
    integer, intent(inout) :: lval(nlev), idx_l_state(nlev)

    integer :: i, j, itmp
    real(dp) :: etmp

    do i = 1, nlev-1
       do j = i+1, nlev
          if (energy(j) < energy(i)) then
             etmp = energy(i)
             energy(i) = energy(j)
             energy(j) = etmp

             itmp = lval(i)
             lval(i) = lval(j)
             lval(j) = itmp

             itmp = idx_l_state(i)
             idx_l_state(i) = idx_l_state(j)
             idx_l_state(j) = itmp
          end if
       end do
    end do
  end subroutine sort_levels


  subroutine sort_levels_with_orbitals(nlev, nr, energy, lval, idx_l_state, u_store)
    integer, intent(in) :: nlev, nr
    real(dp), intent(inout) :: energy(nlev)
    integer, intent(inout) :: lval(nlev), idx_l_state(nlev)
    real(dp), intent(inout) :: u_store(nr,nlev)

    integer :: i, j, itmp, k
    real(dp) :: etmp, utmp

    do i = 1, nlev-1
       do j = i+1, nlev
          if (energy(j) < energy(i)) then
             etmp = energy(i)
             energy(i) = energy(j)
             energy(j) = etmp

             itmp = lval(i)
             lval(i) = lval(j)
             lval(j) = itmp

             itmp = idx_l_state(i)
             idx_l_state(i) = idx_l_state(j)
             idx_l_state(j) = itmp

             do k = 1, nr
                utmp = u_store(k,i)
                u_store(k,i) = u_store(k,j)
                u_store(k,j) = utmp
             end do
          end if
       end do
    end do
  end subroutine sort_levels_with_orbitals

end module sorting_module

module ks_density_module
  use precision_module
  use constants_module
  use radial_grid_module, only: integrate_r3
  implicit none
contains

  subroutine initialize_fermi_density(nr, r, w, A, Nq, radius0, diffuseness, rho_q)
    integer, intent(in) :: nr, A, Nq
    real(dp), intent(in) :: r(nr), w(nr), radius0, diffuseness
    real(dp), intent(out) :: rho_q(nr)

    real(dp) :: R0, a0, norm
    integer :: i

    if (radius0 > 0.0_dp) then
       R0 = radius0
    else
       R0 = 1.20_dp * real(A,dp)**(1.0_dp/3.0_dp)
    end if

    if (diffuseness > 0.0_dp) then
       a0 = diffuseness
    else
       a0 = 0.55_dp
    end if

    do i = 1, nr
       rho_q(i) = 1.0_dp / (1.0_dp + exp((r(i)-R0)/a0))
    end do

    norm = integrate_r3(nr, r, w, rho_q)
    if (norm <= 0.0_dp) stop 'ERROR: initial density normalization failed.'

    rho_q(:) = rho_q(:) * real(Nq, dp) / norm
  end subroutine initialize_fermi_density


  subroutine assign_occupations(nlev, energy, lval, Nq, occ)
    integer, intent(in) :: nlev, lval(nlev), Nq
    real(dp), intent(in) :: energy(nlev)
    real(dp), intent(out) :: occ(nlev)

    integer :: a, remaining, degeneracy, fill

    occ(:) = 0.0_dp
    remaining = Nq

    do a = 1, nlev
       if (remaining <= 0) exit
       degeneracy = 2 * (2*lval(a) + 1)
       fill = min(remaining, degeneracy)
       occ(a) = real(fill, dp)
       remaining = remaining - fill
    end do

    if (remaining > 0) then
       write(*,*) 'ERROR: not enough levels to fill Nq particles. Remaining = ', remaining
       stop
    end if
  end subroutine assign_occupations


  subroutine build_density_from_levels(nr, r, nlev, occ, u_store, rho_q)
    integer, intent(in) :: nr, nlev
    real(dp), intent(in) :: r(nr), occ(nlev), u_store(nr,nlev)
    real(dp), intent(out) :: rho_q(nr)

    integer :: i, a
    real(dp) :: rr

    rho_q(:) = 0.0_dp
    do a = 1, nlev
       if (occ(a) <= 0.0_dp) cycle
       do i = 1, nr
          rr = max(r(i), 1.0e-10_dp)
          rho_q(i) = rho_q(i) + occ(a) * u_store(i,a)**2 / (4.0_dp*pi*rr*rr)
       end do
    end do
  end subroutine build_density_from_levels


  real(dp) function density_residual(nr, r, w, rho_old, rho_new, Nq)
    integer, intent(in) :: nr, Nq
    real(dp), intent(in) :: r(nr), w(nr), rho_old(nr), rho_new(nr)
    real(dp) :: diff2(nr)
    integer :: i

    do i = 1, nr
       diff2(i) = (rho_new(i) - rho_old(i))**2
    end do

    density_residual = sqrt(integrate_r3(nr, r, w, diff2) / max(real(Nq,dp), 1.0_dp))
  end function density_residual

end module ks_density_module

module ks_nlf_module
  use precision_module
  use constants_module
  use radial_grid_module, only: first_derivative_uniform
  implicit none
contains

  subroutine compute_ks_tau_from_levels(nr, r, dr, nlev, occ, lval, u_store, tau_ks)
    integer, intent(in) :: nr, nlev
    real(dp), intent(in) :: r(nr), dr, occ(nlev), u_store(nr,nlev)
    integer, intent(in) :: lval(nlev)
    real(dp), intent(out) :: tau_ks(nr)

    real(dp), allocatable :: Rwave(:), dRwave(:)
    real(dp) :: rr, lfac, pref
    integer :: i, a, l

    allocate(Rwave(nr), dRwave(nr))
    tau_ks(:) = 0.0_dp

    do a = 1, nlev
       if (occ(a) <= 0.0_dp) cycle

       l = lval(a)
       lfac = real(l*(l+1), dp)
       pref = occ(a) / (4.0_dp*pi)

       do i = 1, nr
          rr = max(r(i), 1.0e-10_dp)
          Rwave(i) = u_store(i,a) / rr
       end do

       call first_derivative_uniform(nr, dr, Rwave, dRwave)

       do i = 1, nr
          rr = max(r(i), 1.0e-10_dp)
          tau_ks(i) = tau_ks(i) + pref * &
               ( dRwave(i)**2 + lfac * Rwave(i)**2 / (rr*rr) )
       end do
    end do

    deallocate(Rwave, dRwave)
  end subroutine compute_ks_tau_from_levels


  subroutine compute_nlf_from_tau(nr, tau_method, tau_tf, tau_vw, nlf)
    integer, intent(in) :: nr
    real(dp), intent(in) :: tau_method(nr), tau_tf(nr), tau_vw(nr)
    real(dp), intent(out) :: nlf(nr)

    integer :: i
    real(dp) :: ratio

    do i = 1, nr
       ratio = (tau_method(i) - tau_vw(i)) / max(tau_tf(i), tiny_denom)
       nlf(i) = 1.0_dp / (1.0_dp + ratio*ratio)
    end do
  end subroutine compute_nlf_from_tau

end module ks_nlf_module

module ks_io_module
  use precision_module
  use constants_module
  implicit none
contains

  subroutine write_density_file(filename, nr, r, rho)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nr
    real(dp), intent(in) :: r(nr), rho(nr)

    integer :: i, unit
    unit = 31
    open(unit=unit, file=trim(filename), status='replace', action='write')
    write(unit,'(A)') '# r[fm] rho[fm^-3]'
    do i = 1, nr
       write(unit,'(2ES24.15)') r(i), rho(i)
    end do
    close(unit)
  end subroutine write_density_file


  subroutine write_levels_file(filename, nlev, energy, lval, occ)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nlev, lval(nlev)
    real(dp), intent(in) :: energy(nlev), occ(nlev)

    integer :: a, unit
    unit = 32
    open(unit=unit, file=trim(filename), status='replace', action='write')
    write(unit,'(A)') '# index l energy[MeV] degeneracy occupation'
    do a = 1, nlev
       write(unit,'(2I8,ES24.15,I8,ES24.15)') a, lval(a), energy(a), &
            2*(2*lval(a)+1), occ(a)
    end do
    close(unit)
  end subroutine write_levels_file


  subroutine write_scf_history_header(unit, filename)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: filename

    open(unit=unit, file=trim(filename), status='replace', action='write')
    write(unit,'(A)') '# iter residual particle_number_species min_level_energy[MeV] max_occ_energy[MeV]'
  end subroutine write_scf_history_header


  subroutine append_scf_history(unit, iter, residual, number_q, e_min, e_fermi)
    integer, intent(in) :: unit, iter
    real(dp), intent(in) :: residual, number_q, e_min, e_fermi

    write(unit,'(I8,4ES24.15)') iter, residual, number_q, e_min, e_fermi
  end subroutine append_scf_history


  subroutine write_ks_nlf_file(filename, nr, r, rho, tau_tf, tau_vw, tau_ks, nlf_ks)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nr
    real(dp), intent(in) :: r(nr), rho(nr), tau_tf(nr), tau_vw(nr), tau_ks(nr), nlf_ks(nr)

    integer :: i, unit

    unit = 34
    open(unit=unit, file=trim(filename), status='replace', action='write')
    write(unit,'(A)') '# r[fm] rho_species[fm^-3] tau_TF tau_vW tau_KS NLF_KS'
    do i = 1, nr
       write(unit,'(6ES24.15)') r(i), rho(i), tau_tf(i), tau_vw(i), tau_ks(i), nlf_ks(i)
    end do
    close(unit)
  end subroutine write_ks_nlf_file

end module ks_io_module

module ks_scf_module
  use precision_module
  use constants_module
  use radial_grid_module
  use skyrme_skp_module
  use ks_orbital_module
  use sorting_module
  use ks_density_module
  use ks_nlf_module
  use ks_io_module
  implicit none
contains

  subroutine run_spherical_ks_skp(A, Nq, nr, rbox, lmax, nstates_per_l, &
       max_scf, scf_tol, mix_alpha, init_radius, init_diffuseness, &
       density_species_file, density_total_file, levels_file, history_file, ks_nlf_file, &
       r, dr, w, rho_q, rho_tot)

    integer, intent(in) :: A, Nq, nr, lmax, nstates_per_l, max_scf
    real(dp), intent(in) :: rbox, scf_tol, mix_alpha, init_radius, init_diffuseness
    character(len=*), intent(in) :: density_species_file, density_total_file
    character(len=*), intent(in) :: levels_file, history_file, ks_nlf_file

    real(dp), allocatable, intent(out) :: r(:), w(:), rho_q(:), rho_tot(:)
    real(dp), intent(out) :: dr

    real(dp), allocatable :: vks(:), rho_new(:), rho_old(:)
    real(dp), allocatable :: eigval_l(:), u_l(:,:)
    real(dp), allocatable :: energy(:), occ(:), u_store(:,:)
    real(dp), allocatable :: tau_tf_ks(:), tau_vw_ks(:), tau_ks(:), nlf_ks(:)
    integer, allocatable :: lval(:), idx_l_state(:)
    integer :: nlev, lev, l, ist, iter, hist_unit
    real(dp) :: res, number_q, e_min, e_fermi

    allocate(r(nr), w(nr), rho_q(nr), rho_tot(nr))
    allocate(vks(nr), rho_new(nr), rho_old(nr))

    call make_uniform_interior_grid(nr, rbox, r, dr, w)
    call initialize_fermi_density(nr, r, w, A, Nq, init_radius, init_diffuseness, rho_q)

    nlev = (lmax + 1) * nstates_per_l
    allocate(energy(nlev), occ(nlev), lval(nlev), idx_l_state(nlev))
    allocate(u_store(nr,nlev))
    allocate(eigval_l(nr), u_l(nr,nr))

    hist_unit = 77
    call write_scf_history_header(hist_unit, history_file)

    write(*,*) '-----------------------------------------------'
    write(*,*) 'Starting spherical KS-SkP SCF'
    write(*,*) 'A, Nq              = ', A, Nq
    write(*,*) 'nr, rbox, dr       = ', nr, rbox, dr
    write(*,*) 'lmax, states/l     = ', lmax, nstates_per_l
    write(*,*) 'max_scf, tol, mix  = ', max_scf, scf_tol, mix_alpha
    write(*,*) '-----------------------------------------------'

    do iter = 1, max_scf
       rho_old(:) = rho_q(:)
       rho_tot(:) = 2.0_dp * rho_q(:)

       call compute_skp_potential(nr, r, dr, rho_tot, vks)

       lev = 0
       do l = 0, lmax
          call solve_radial_l(nr, r, dr, l, vks, eigval_l, u_l)
          do ist = 1, nstates_per_l
             lev = lev + 1
             energy(lev)      = eigval_l(ist)
             lval(lev)        = l
             idx_l_state(lev) = ist
             u_store(:,lev)   = u_l(:,ist)
          end do
       end do

       call sort_levels_with_orbitals(nlev, nr, energy, lval, idx_l_state, u_store)

       call assign_occupations(nlev, energy, lval, Nq, occ)
       call build_density_from_levels(nr, r, nlev, occ, u_store, rho_new)

       res = density_residual(nr, r, w, rho_old, rho_new, Nq)

       rho_q(:) = (1.0_dp - mix_alpha) * rho_old(:) + mix_alpha * rho_new(:)
       where (rho_q < 0.0_dp) rho_q = 0.0_dp

       number_q = integrate_r3(nr, r, w, rho_q)

       ! Renormalize after mixing to exactly Nq particles.
       if (number_q > 0.0_dp) then
          rho_q(:) = rho_q(:) * real(Nq,dp) / number_q
       end if
       number_q = integrate_r3(nr, r, w, rho_q)

       e_min = energy(1)
       e_fermi = highest_occupied_energy(nlev, energy, occ)

       call append_scf_history(hist_unit, iter, res, number_q, e_min, e_fermi)

       if (mod(iter, 10) == 0 .or. iter == 1 .or. res < scf_tol) then
          write(*,'(A,I6,A,ES12.4,A,F12.6,A,F12.6)') &
               'KS iter ', iter, ' residual=', res, ' Nq=', number_q, ' eF=', e_fermi
       end if

       if (res < scf_tol) exit
    end do

    close(hist_unit)

    rho_tot(:) = 2.0_dp * rho_q(:)

    ! Recompute the final KS spectrum/orbitals from the final mixed density.
    ! This makes the output levels and the KS-NLF black curve consistent with rho_q.
    call compute_skp_potential(nr, r, dr, rho_tot, vks)

    lev = 0
    do l = 0, lmax
       call solve_radial_l(nr, r, dr, l, vks, eigval_l, u_l)
       do ist = 1, nstates_per_l
          lev = lev + 1
          energy(lev)      = eigval_l(ist)
          lval(lev)        = l
          idx_l_state(lev) = ist
          u_store(:,lev)   = u_l(:,ist)
       end do
    end do

    call sort_levels_with_orbitals(nlev, nr, energy, lval, idx_l_state, u_store)
    call assign_occupations(nlev, energy, lval, Nq, occ)

    allocate(tau_tf_ks(nr), tau_vw_ks(nr), tau_ks(nr), nlf_ks(nr))
    call compute_tf_vw_tau(nr, r, dr, rho_q, tau_tf_ks, tau_vw_ks)
    call compute_ks_tau_from_levels(nr, r, dr, nlev, occ, lval, u_store, tau_ks)
    call compute_nlf_from_tau(nr, tau_ks, tau_tf_ks, tau_vw_ks, nlf_ks)
    call write_ks_nlf_file(ks_nlf_file, nr, r, rho_q, tau_tf_ks, tau_vw_ks, tau_ks, nlf_ks)

    write(*,*) 'KS NLF written to         : ', trim(ks_nlf_file)

    deallocate(tau_tf_ks, tau_vw_ks, tau_ks, nlf_ks)

    call write_density_file(density_species_file, nr, r, rho_q)
    call write_density_file(density_total_file,   nr, r, rho_tot)
    call write_levels_file(levels_file, nlev, energy, lval, occ)

    write(*,*) 'KS-SCF finished.'
    write(*,*) 'Species density written to: ', trim(density_species_file)
    write(*,*) 'Total density written to  : ', trim(density_total_file)
    write(*,*) 'Levels written to         : ', trim(levels_file)

    deallocate(vks, rho_new, rho_old, eigval_l, u_l, energy, occ, lval, idx_l_state, u_store)
  end subroutine run_spherical_ks_skp


  subroutine rebuild_sorted_orbitals(nr, r, dr, lmax, nstates_per_l, vks, &
       nlev, energy_sorted, lval_sorted, idx_l_state_sorted, u_store_sorted)
    integer, intent(in) :: nr, lmax, nstates_per_l, nlev
    real(dp), intent(in) :: r(nr), dr, vks(nr), energy_sorted(nlev)
    integer, intent(in) :: lval_sorted(nlev), idx_l_state_sorted(nlev)
    real(dp), intent(out) :: u_store_sorted(nr,nlev)

    real(dp), allocatable :: eigval_l(:), u_l(:,:)
    integer :: l, a, ist

    allocate(eigval_l(nr), u_l(nr,nr))

    do a = 1, nlev
       l = lval_sorted(a)
       ist = idx_l_state_sorted(a)
       call solve_radial_l(nr, r, dr, l, vks, eigval_l, u_l)
       u_store_sorted(:,a) = u_l(:,ist)
       ! energy_sorted(a) should match eigval_l(ist), but we do not enforce it here.
    end do

    deallocate(eigval_l, u_l)
  end subroutine rebuild_sorted_orbitals


  real(dp) function highest_occupied_energy(nlev, energy, occ)
    integer, intent(in) :: nlev
    real(dp), intent(in) :: energy(nlev), occ(nlev)
    integer :: a

    highest_occupied_energy = energy(1)
    do a = 1, nlev
       if (occ(a) > 0.0_dp) highest_occupied_energy = energy(a)
    end do
  end function highest_occupied_energy

end module ks_scf_module

module density_file_module
  use precision_module
  implicit none
contains

  subroutine read_density_2col(filename, density_scale, n, r, rho)
    character(len=*), intent(in) :: filename
    real(dp), intent(in) :: density_scale
    integer, intent(out) :: n
    real(dp), allocatable, intent(out) :: r(:), rho(:)

    integer :: unit, ios, count, i
    character(len=512) :: line, sline
    real(dp) :: x, y

    unit = 41
    count = 0

    open(unit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(*,*) 'ERROR: cannot open density file: ', trim(filename)
       stop
    end if

    do
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) exit
       sline = adjustl(line)
       if (len_trim(sline) == 0) cycle
       if (sline(1:1) == '#') cycle
       if (sline(1:1) == '!') cycle
       read(sline, *, iostat=ios) x, y
       if (ios == 0) count = count + 1
    end do
    close(unit)

    if (count < 3) stop 'ERROR: too few density data points.'

    n = count
    allocate(r(n), rho(n))

    open(unit=unit, file=trim(filename), status='old', action='read')
    i = 0
    do
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) exit
       sline = adjustl(line)
       if (len_trim(sline) == 0) cycle
       if (sline(1:1) == '#') cycle
       if (sline(1:1) == '!') cycle
       read(sline, *, iostat=ios) x, y
       if (ios == 0) then
          i = i + 1
          r(i)   = x
          rho(i) = max(density_scale*y, 0.0_dp)
       end if
    end do
    close(unit)
  end subroutine read_density_2col

end module density_file_module

module gauss_module
  use precision_module
  use constants_module, only: pi
  implicit none
contains

  subroutine gauss_legendre(x1, x2, x, w, n)
    real(dp), intent(in) :: x1, x2
    integer, intent(in) :: n
    real(dp), intent(out) :: x(n), w(n)

    integer :: m, j, i
    real(dp) :: xm, xl, z, z1, p1, p2, p3, pp
    real(dp), parameter :: eps = 1.0e-14_dp

    m = (n + 1)/2
    xm = 0.5_dp * (x2 + x1)
    xl = 0.5_dp * (x2 - x1)

    do i = 1, m
       z = cos(pi * (real(i,dp) - 0.25_dp) / (real(n,dp) + 0.5_dp))

       do
          p1 = 1.0_dp
          p2 = 0.0_dp
          do j = 1, n
             p3 = p2
             p2 = p1
             p1 = ((2.0_dp*real(j,dp)-1.0_dp)*z*p2 - &
                   (real(j,dp)-1.0_dp)*p3) / real(j,dp)
          end do
          pp = real(n,dp) * (z*p1 - p2) / (z*z - 1.0_dp)
          z1 = z
          z = z1 - p1/pp
          if (abs(z-z1) < eps) exit
       end do

       x(i)       = xm - xl*z
       x(n+1-i)   = xm + xl*z
       w(i)       = 2.0_dp*xl / ((1.0_dp - z*z)*pp*pp)
       w(n+1-i)   = w(i)
    end do
  end subroutine gauss_legendre

end module gauss_module

module lindhard_kernel_module
  use precision_module
  use constants_module
  implicit none

  ! Dimensionless kernel table for the density-dependent-kF reference calculation:
  ! f1(R; kF) = kF^3 * g(kF R),  g(y) = f1_of_R(R=y, kf=1, qmax=xqmax_dimless).
  logical :: g_ready = .false.
  real(dp), allocatable :: g_m2tab(:)   ! natural-spline second derivatives of g

  ! Fitted nonlocal correction kernel (added onto the constant-kF f1
  ! table; two-column file R[fm], delta_f).  Because the correction enters as a
  ! pure symmetric two-point kernel plus a c_CORR shift, the exact variational
  ! machinery (derivatives, checker, ITE, minimizer) is inherited without new
  ! derivative code.
  logical :: extra_kernel_loaded = .false.
  integer :: nek = 0
  real(dp), allocatable :: ek_R(:), ek_f(:)

  ! Fitted correction to the DIMENSIONLESS dd kernel shape g(y).  The table
  ! (y, delta_g) is added onto g_gtab inside build_dimensionless_f1_table
  ! BEFORE the spline is built, so the exact variational machinery
  ! (df1/dt = g + (y/3) g') carries over without new derivative code.
  logical :: extra_g_loaded = .false.
  integer :: neg = 0
  real(dp), allocatable :: eg_y(:), eg_g(:)
  real(dp) :: g_ymax_saved = -1.0_dp
  ! Full cache key for the dimensionless table: a stale table
  ! must not be reused when xqmax, nq, or the loaded delta_g change even if
  ! (ny, ymax) are unchanged.
  real(dp) :: g_xqmax_saved = -1.0_dp
  integer  :: g_nq_saved = -1
  real(dp) :: g_dghash_saved = -9.99e30_dp
  real(dp), allocatable :: g_ytab(:), g_gtab(:)

contains

  real(dp) function lindhard_l(x)
    real(dp), intent(in) :: x
    real(dp) :: ax, logterm

    ax = abs(x)
    if (ax < 1.0e-8_dp) then
       lindhard_l = 1.0_dp - x*x/3.0_dp - x**4/15.0_dp
    else if (abs(ax - 1.0_dp) < 1.0e-7_dp) then
       lindhard_l = 0.5_dp
    else
       logterm = log(abs((1.0_dp + ax)/(1.0_dp - ax)))
       lindhard_l = 0.5_dp * &
            (1.0_dp + (1.0_dp - ax*ax)/(2.0_dp*ax) * logterm)
    end if
  end function lindhard_l


  real(dp) function lindhard_l1(x)
    real(dp), intent(in) :: x
    real(dp) :: lx

    lx = lindhard_l(x)
    lindhard_l1 = (5.0_dp/8.0_dp) * &
         (1.0_dp/lx - 3.0_dp*x*x + 3.0_dp/5.0_dp)
  end function lindhard_l1


  real(dp) function f1_of_R(R, kf, qmax, nq)
    real(dp), intent(in) :: R, kf, qmax
    integer, intent(in) :: nq

    integer :: iq
    real(dp) :: hq, q, x, val, weight, sum

    hq = qmax / real(nq, dp)
    sum = 0.0_dp

    do iq = 0, nq
       q = real(iq, dp) * hq
       x = q / (2.0_dp * kf)

       if (iq == 0 .or. iq == nq) then
          weight = 1.0_dp
       else if (mod(iq, 2) == 1) then
          weight = 4.0_dp
       else
          weight = 2.0_dp
       end if

       if (abs(R) < 1.0e-10_dp) then
          val = q*q * lindhard_l1(x)
       else
          val = q * lindhard_l1(x) * sin(q*R)
       end if

       sum = sum + weight * val
    end do

    sum = sum * hq / 3.0_dp

    if (abs(R) < 1.0e-10_dp) then
       f1_of_R = sum / (2.0_dp*pi*pi)
    else
       f1_of_R = sum / (2.0_dp*pi*pi*R)
    end if
  end function f1_of_R


  subroutine load_extra_kernel(fname)
    character(len=*), intent(in) :: fname
    integer :: unit, ios, n
    real(dp) :: a, b
    character(len=256) :: line
    open(newunit=unit, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'ERROR: cannot open extra_kernel_file.'
    n = 0
    do
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) exit
       if (len_trim(line) == 0) cycle
       if (line(1:1) == '#') cycle
       n = n + 1
    end do
    rewind(unit)
    if (allocated(ek_R)) deallocate(ek_R, ek_f)
    allocate(ek_R(n), ek_f(n))
    nek = 0
    do
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) exit
       if (len_trim(line) == 0) cycle
       if (line(1:1) == '#') cycle
       read(line, *) a, b
       nek = nek + 1
       ek_R(nek) = a
       ek_f(nek) = b
    end do
    close(unit)
    extra_kernel_loaded = .true.
    write(*,'(A,I6,A,A)') ' Extra fitted kernel loaded: ', nek, ' points from ', trim(fname)
  end subroutine load_extra_kernel


  real(dp) function interp_extra_kernel(R)
    real(dp), intent(in) :: R
    integer :: i
    real(dp) :: t
    if (R <= ek_R(1)) then
       interp_extra_kernel = ek_f(1)
       return
    end if
    if (R >= ek_R(nek)) then
       interp_extra_kernel = 0.0_dp
       return
    end if
    i = int((R - ek_R(1)) / (ek_R(2) - ek_R(1))) + 1
    if (i < 1) i = 1
    if (i >= nek) i = nek - 1
    t = (R - ek_R(i)) / (ek_R(i+1) - ek_R(i))
    interp_extra_kernel = (1.0_dp - t)*ek_f(i) + t*ek_f(i+1)
  end function interp_extra_kernel


  subroutine load_extra_g(fname)
    character(len=*), intent(in) :: fname
    integer :: unit, ios, n
    real(dp) :: a, b
    character(len=256) :: line
    open(newunit=unit, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'ERROR: cannot open extra_g_file.'
    n = 0
    do
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) exit
       if (len_trim(line) == 0) cycle
       if (line(1:1) == '#') cycle
       n = n + 1
    end do
    rewind(unit)
    if (allocated(eg_y)) deallocate(eg_y, eg_g)
    allocate(eg_y(n), eg_g(n))
    neg = 0
    do
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) exit
       if (len_trim(line) == 0) cycle
       if (line(1:1) == '#') cycle
       read(line, *) a, b
       neg = neg + 1
       eg_y(neg) = a
       eg_g(neg) = b
    end do
    close(unit)
    extra_g_loaded = .true.
    write(*,'(A,I6,A,A)') ' Extra fitted delta-g loaded: ', neg, ' points from ', trim(fname)
  end subroutine load_extra_g


  real(dp) function interp_extra_g(y)
    real(dp), intent(in) :: y
    integer :: i
    real(dp) :: tt
    if (y <= eg_y(1)) then
       interp_extra_g = eg_g(1)
       return
    end if
    if (y >= eg_y(neg)) then
       interp_extra_g = 0.0_dp
       return
    end if
    i = int((y - eg_y(1)) / (eg_y(2) - eg_y(1))) + 1
    if (i < 1) i = 1
    if (i >= neg) i = neg - 1
    tt = (y - eg_y(i)) / (eg_y(i+1) - eg_y(i))
    interp_extra_g = (1.0_dp - tt)*eg_g(i) + tt*eg_g(i+1)
  end function interp_extra_g


  subroutine build_f1_table(kf, qmax, nq_input, nr_tab, rmax_tab, rtab, ftab)
    real(dp), intent(in) :: kf, qmax, rmax_tab
    integer, intent(in) :: nq_input, nr_tab
    real(dp), intent(out) :: rtab(nr_tab), ftab(nr_tab)

    integer :: i, nq
    real(dp) :: dR

    nq = nq_input
    if (mod(nq,2) /= 0) nq = nq + 1

    dR = rmax_tab / real(nr_tab-1, dp)
    do i = 1, nr_tab
       rtab(i) = real(i-1, dp) * dR
       ftab(i) = f1_of_R(rtab(i), kf, qmax, nq)
    end do
    if (extra_kernel_loaded) then
       block
         integer :: iek
         do iek = 1, nr_tab
            ftab(iek) = ftab(iek) + interp_extra_kernel(rtab(iek))
         end do
       end block
    end if

  end subroutine build_f1_table


  real(dp) function interp_f1(R, nr_tab, rtab, ftab)
    real(dp), intent(in) :: R
    integer, intent(in) :: nr_tab
    real(dp), intent(in) :: rtab(nr_tab), ftab(nr_tab)

    integer :: i
    real(dp) :: dR, t

    if (R <= rtab(1)) then
       interp_f1 = ftab(1)
       return
    end if
    if (R >= rtab(nr_tab)) then
       interp_f1 = 0.0_dp
       return
    end if

    dR = rtab(2) - rtab(1)
    i = int((R - rtab(1))/dR) + 1
    if (i < 1) i = 1
    if (i >= nr_tab) i = nr_tab - 1

    t = (R - rtab(i)) / (rtab(i+1) - rtab(i))
    interp_f1 = (1.0_dp - t)*ftab(i) + t*ftab(i+1)
  end function interp_f1


  subroutine check_f1_norm(nr_tab, rtab, ftab)
    integer, intent(in) :: nr_tab
    real(dp), intent(in) :: rtab(nr_tab), ftab(nr_tab)

    real(dp) :: norm, dR, w
    integer :: i

    dR = rtab(2) - rtab(1)
    norm = 0.0_dp
    do i = 1, nr_tab
       if (i == 1 .or. i == nr_tab) then
          w = 0.5_dp*dR
       else
          w = dR
       end if
       norm = norm + 4.0_dp*pi*rtab(i)**2 * ftab(i) * w
    end do

    write(*,'(A,ES18.8)') 'Check int 4*pi*R^2*f1(R)dR = ', norm
  end subroutine check_f1_norm


  subroutine build_dimensionless_f1_table(xqmax, nq_input, ny, ymax)
    ! Tabulate g(y) = f1_of_R(y, kf=1, xqmax) on y in [0, ymax], so that
    ! f1(R; kF) = kF^3 * g(kF R) can be evaluated for arbitrary kF.
    real(dp), intent(in) :: xqmax, ymax
    integer, intent(in) :: nq_input, ny

    integer :: i, nq
    real(dp) :: dy, norm, ww
    real(dp) :: dghash

    nq = nq_input
    if (mod(nq,2) /= 0) nq = nq + 1

    ! content hash of the loaded delta_g table (0 if none), so a change in the
    ! fitted correction invalidates a cached kernel table.
    dghash = 0.0_dp
    if (extra_g_loaded) then
       do i = 1, size(eg_g)
          dghash = dghash + real(i,dp)*eg_g(i) + eg_y(i)
       end do
    end if

    if (g_ready .and. allocated(g_ytab)) then
       if (size(g_ytab) == ny .and. abs(g_ymax_saved - ymax) < 1.0e-12_dp &
           .and. abs(g_xqmax_saved - xqmax) < 1.0e-12_dp .and. g_nq_saved == nq &
           .and. abs(g_dghash_saved - dghash) < 1.0e-12_dp) return
    end if

    if (allocated(g_ytab)) deallocate(g_ytab)
    if (allocated(g_gtab)) deallocate(g_gtab)
    allocate(g_ytab(ny), g_gtab(ny))

    dy = ymax / real(ny-1, dp)
    do i = 1, ny
       g_ytab(i) = real(i-1, dp) * dy
       g_gtab(i) = f1_of_R(g_ytab(i), 1.0_dp, xqmax, nq)
    end do
    if (extra_g_loaded) then
       do i = 1, ny
          g_gtab(i) = g_gtab(i) + interp_extra_g(g_ytab(i))
       end do
       write(*,'(A)') ' delta-g added to the dimensionless dd kernel table.'
    end if
    g_ymax_saved = ymax
    g_xqmax_saved = xqmax
    g_nq_saved = nq
    g_dghash_saved = dghash
    call build_g_spline()
    g_ready = .true.

    norm = 0.0_dp
    do i = 1, ny
       if (i == 1 .or. i == ny) then
          ww = 0.5_dp*dy
       else
          ww = dy
       end if
       norm = norm + 4.0_dp*pi*g_ytab(i)**2 * g_gtab(i) * ww
    end do
    write(*,'(A,ES18.8)') 'Dimensionless kernel check int 4*pi*y^2*g(y)dy = ', norm
  end subroutine build_dimensionless_f1_table


  subroutine build_g_spline()
    ! Natural cubic spline second derivatives M_i of g on the uniform y grid,
    ! solved by the Thomas algorithm.  A C^2 interpolant makes the discrete
    ! dd-kernel energy a C^2 function of the density, which is required for the
    ! finite-difference derivative checker to exhibit clean O(eps^2) scaling
    ! (a piecewise-LINEAR interpolant has slope kinks at every knot; with
    ! ~nr^2*nmu pair distances straddling knots the FD error is dominated by
    ! kink crossings and never reaches roundoff).
    integer :: i, ny
    real(dp) :: h
    real(dp), allocatable :: a(:), b(:), c(:), d(:)

    ny = size(g_ytab)
    if (allocated(g_m2tab)) deallocate(g_m2tab)
    allocate(g_m2tab(ny), a(ny), b(ny), c(ny), d(ny))
    h = g_ytab(2) - g_ytab(1)

    ! natural BCs: M_1 = M_ny = 0
    b(1) = 1.0_dp; c(1) = 0.0_dp; d(1) = 0.0_dp
    do i = 2, ny-1
       a(i) = h/6.0_dp
       b(i) = 2.0_dp*h/3.0_dp
       c(i) = h/6.0_dp
       d(i) = (g_gtab(i+1) - 2.0_dp*g_gtab(i) + g_gtab(i-1)) / h
    end do
    a(ny) = 0.0_dp; b(ny) = 1.0_dp; d(ny) = 0.0_dp

    ! Thomas sweep
    do i = 2, ny
       d(i) = d(i) - (a(i)/b(i-1))*d(i-1)
       b(i) = b(i) - (a(i)/b(i-1))*c(i-1)
    end do
    g_m2tab(ny) = d(ny)/b(ny)
    do i = ny-1, 1, -1
       g_m2tab(i) = (d(i) - c(i)*g_m2tab(i+1)) / b(i)
    end do
    deallocate(a, b, c, d)
  end subroutine build_g_spline


  real(dp) function interp_g_dimless(y)
    real(dp), intent(in) :: y
    real(dp) :: gv, gs
    call interp_g_slope(y, gv, gs)
    interp_g_dimless = gv
  end function interp_g_dimless


  subroutine interp_g_slope(y, gval, gslope)
    ! Piecewise-linear interpolant of g(y) together with the exact slope of the
    ! interpolant on the current segment.  Using the segment slope (rather than a
    ! smooth reconstruction) keeps  d f1 / d t  the exact derivative of the
    ! discrete energy that actually uses this interpolant, so the finite
    ! difference checker can reach roundoff level in dd-kernel mode as well.
    real(dp), intent(in) :: y
    real(dp), intent(out) :: gval, gslope
    integer :: i, ny
    real(dp) :: dy, t

    if (.not. g_ready) stop 'ERROR: dimensionless f1 table is not built.'
    ny = size(g_ytab)
    if (y <= g_ytab(1)) then
       gval = g_gtab(1); gslope = 0.0_dp
       return
    end if
    if (y >= g_ytab(ny)) then
       gval = 0.0_dp; gslope = 0.0_dp
       return
    end if
    dy = g_ytab(2) - g_ytab(1)
    i = int((y - g_ytab(1))/dy) + 1
    if (i < 1) i = 1
    if (i >= ny) i = ny - 1
    block
      real(dp) :: A, B
      A = (g_ytab(i+1) - y)/dy
      B = 1.0_dp - A
      gval = A*g_gtab(i) + B*g_gtab(i+1) &
           + ((A**3 - A)*g_m2tab(i) + (B**3 - B)*g_m2tab(i+1))*dy*dy/6.0_dp
      gslope = (g_gtab(i+1) - g_gtab(i))/dy &
           - (3.0_dp*A*A - 1.0_dp)*dy/6.0_dp*g_m2tab(i) &
           + (3.0_dp*B*B - 1.0_dp)*dy/6.0_dp*g_m2tab(i+1)
    end block
  end subroutine interp_g_slope

end module lindhard_kernel_module

module nlf_postprocess_module
  use precision_module
  use constants_module
  use radial_grid_module, only: compute_tf_vw_tau
  use gauss_module
  use lindhard_kernel_module
  implicit none

  ! Cached constant-kF convolution matrix:
  !   conv(i) = sum_j Mconv(i,j) * rho(j)^(5/6),
  !   Mconv(i,j) = 2 pi r_j^2 w_j sum_k wmu_k f1(|r_i - r_j|_k),
  ! rebuilt automatically whenever the grid, nmu, or the kernel table changes.
  logical :: mconv_ready = .false.
  integer :: mconv_nr = -1, mconv_nmu = -1
  real(dp) :: mconv_sig = 0.0_dp
  real(dp), allocatable :: Mconv(:,:)

contains

  subroutine compute_constkf_convolution(nr, r, w, rho, nmu, nr_tab, rtab, ftab, conv)
    integer, intent(in) :: nr, nmu, nr_tab
    real(dp), intent(in) :: r(nr), w(nr), rho(nr)
    real(dp), intent(in) :: rtab(nr_tab), ftab(nr_tab)
    real(dp), intent(out) :: conv(nr)

    real(dp), allocatable :: g(:)
    integer :: i, j
    real(dp) :: sig, acc

    ! signature of (grid, quadrature, kernel table)
    sig = real(nr,dp) + real(nmu,dp)*1.0e3_dp + r(nr) + w(1) &
        + sum(ftab)*1.0e-3_dp + ftab(1) + rtab(nr_tab)
    if (.not. (mconv_ready .and. mconv_nr == nr .and. mconv_nmu == nmu .and. &
         abs(mconv_sig - sig) < 1.0e-10_dp*max(1.0_dp, abs(sig)))) then
       call build_constkf_convolution_matrix(nr, r, w, nmu, nr_tab, rtab, ftab)
       mconv_nr = nr
       mconv_nmu = nmu
       mconv_sig = sig
       mconv_ready = .true.
    end if

    allocate(g(nr))
    do j = 1, nr
       g(j) = max(rho(j), 0.0_dp)**(5.0_dp/6.0_dp)
    end do

    do i = 1, nr
       acc = 0.0_dp
       do j = 1, nr
          acc = acc + Mconv(i,j) * g(j)
       end do
       conv(i) = acc
    end do
    deallocate(g)
  end subroutine compute_constkf_convolution


  subroutine build_constkf_convolution_matrix(nr, r, w, nmu, nr_tab, rtab, ftab)
    integer, intent(in) :: nr, nmu, nr_tab
    real(dp), intent(in) :: r(nr), w(nr)
    real(dp), intent(in) :: rtab(nr_tab), ftab(nr_tab)

    real(dp), allocatable :: mu(:), wmu(:)
    integer :: i, j, k
    real(dp) :: ri, rj, Rij, angular_sum, kval

    write(*,*) 'Building constant-kF convolution matrix cache ...'
    if (allocated(Mconv)) deallocate(Mconv)
    allocate(Mconv(nr,nr))
    allocate(mu(nmu), wmu(nmu))
    call gauss_legendre(-1.0_dp, 1.0_dp, mu, wmu, nmu)

    do i = 1, nr
       ri = r(i)
       do j = 1, nr
          rj = r(j)
          angular_sum = 0.0_dp
          do k = 1, nmu
             Rij = sqrt(max(0.0_dp, ri*ri + rj*rj - 2.0_dp*ri*rj*mu(k)))
             kval = interp_f1(Rij, nr_tab, rtab, ftab)
             angular_sum = angular_sum + wmu(k) * kval
          end do
          Mconv(i,j) = 2.0_dp*pi * rj*rj * w(j) * angular_sum
       end do
    end do
    deallocate(mu, wmu)
  end subroutine build_constkf_convolution_matrix


  subroutine compute_tau_nl_and_nlf(n, rho, conv, tau_tf, tau_vw, tau_nl, nlf)
    integer, intent(in) :: n
    real(dp), intent(in) :: rho(n), conv(n), tau_tf(n), tau_vw(n)
    real(dp), intent(out) :: tau_nl(n), nlf(n)

    integer :: i
    real(dp) :: rh, ratio

    do i = 1, n
       rh = max(rho(i), rho_floor)
       tau_nl(i) = tau_vw(i) &
            + c_NL   * rh**(5.0_dp/6.0_dp) * conv(i) &
            - c_CORR * rh**(5.0_dp/3.0_dp)

       ratio = (tau_nl(i) - tau_vw(i)) / max(tau_tf(i), tiny_denom)
       nlf(i) = 1.0_dp / (1.0_dp + ratio*ratio)
    end do
  end subroutine compute_tau_nl_and_nlf


  subroutine run_nlf_constkf(nr, r, dr, w, rho, kf_const, qmax, nq, nmu, nr_tab, &
       nlf_file, kernel_file)
    integer, intent(in) :: nr, nq, nmu, nr_tab
    real(dp), intent(in) :: r(nr), dr, w(nr), rho(nr), kf_const, qmax
    character(len=*), intent(in) :: nlf_file, kernel_file

    real(dp), allocatable :: rtab(:), ftab(:)
    real(dp), allocatable :: conv(:), tau_tf(:), tau_vw(:), tau_nl(:), nlf(:)
    real(dp) :: rtab_max

    allocate(rtab(nr_tab), ftab(nr_tab))
    allocate(conv(nr), tau_tf(nr), tau_vw(nr), tau_nl(nr), nlf(nr))

    rtab_max = 2.0_dp * r(nr)

    write(*,*) '-----------------------------------------------'
    write(*,*) 'Starting NL const-kF NLF postprocessing'
    write(*,*) 'kf_const [fm^-1] = ', kf_const
    write(*,*) 'qmax [fm^-1]     = ', qmax
    write(*,*) 'nq, nmu, nr_tab  = ', nq, nmu, nr_tab
    write(*,*) '-----------------------------------------------'

    call build_f1_table(kf_const, qmax, nq, nr_tab, rtab_max, rtab, ftab)
    call check_f1_norm(nr_tab, rtab, ftab)
    call write_kernel_file(kernel_file, nr_tab, rtab, ftab)

    call compute_tf_vw_tau(nr, r, dr, rho, tau_tf, tau_vw)
    call compute_constkf_convolution(nr, r, w, rho, nmu, nr_tab, rtab, ftab, conv)
    call compute_tau_nl_and_nlf(nr, rho, conv, tau_tf, tau_vw, tau_nl, nlf)
    call write_nlf_file(nlf_file, nr, r, rho, conv, tau_tf, tau_vw, tau_nl, nlf)

    write(*,*) 'NLF output written to   : ', trim(nlf_file)
    write(*,*) 'Kernel output written to: ', trim(kernel_file)

    deallocate(rtab, ftab, conv, tau_tf, tau_vw, tau_nl, nlf)
  end subroutine run_nlf_constkf



  subroutine write_nlf_file(filename, n, r, rho, conv, tau_tf, tau_vw, tau_nl, nlf)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: n
    real(dp), intent(in) :: r(n), rho(n), conv(n), tau_tf(n), tau_vw(n), tau_nl(n), nlf(n)

    integer :: i, unit
    unit = 51
    open(unit=unit, file=trim(filename), status='replace', action='write')
    write(unit,'(A)') '# r[fm] rho_species[fm^-3] conv tau_TF tau_vW tau_NL NLF_NL'
    do i = 1, n
       write(unit,'(7ES24.15)') r(i), rho(i), conv(i), tau_tf(i), tau_vw(i), tau_nl(i), nlf(i)
    end do
    close(unit)
  end subroutine write_nlf_file


  subroutine write_kernel_file(filename, n, rtab, ftab)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: n
    real(dp), intent(in) :: rtab(n), ftab(n)

    integer :: i, unit
    unit = 52
    open(unit=unit, file=trim(filename), status='replace', action='write')
    write(unit,'(A)') '# R[fm] f1_const_kF[R] [fm^-3]'
    do i = 1, n
       write(unit,'(2ES24.15)') rtab(i), ftab(i)
    end do
    close(unit)
  end subroutine write_kernel_file

end module nlf_postprocess_module

module semilocal_correction_module
  use precision_module
  use constants_module
  use radial_grid_module, only: first_derivative_uniform, second_derivative_uniform, integrate_r3
  implicit none

  ! Fitted semilocal correction functional:
  !
  !   T_SL[rho] = C int rho^(5/3) S(rho) F(p_eff, q_eff) d^3r ,   C = hbar^2/(2m)
  !
  ! with the semilocal invariants
  !   p = |grad rho|^2 / (a rho^(8/3)),  q = lap rho / (a rho^(5/3)),  a = 4 (3 pi^2)^(2/3),
  ! smooth clipping that bounds the dimensionless invariants,
  !   p_eff = pclip * tanh(p/pclip),  q_eff = qclip * tanh(q/qclip),
  ! a smooth low-density switch that turns the correction off in the far tail
  !   S(rho) = rho^2 / (rho^2 + rho_switch^2),
  ! and a polynomial basis obeying F(0,0)=0 so uniform matter is unchanged:
  !   F = c10 p + c20 p^2 + c11 p q + c02 q^2 + c30 p^3 + c21 p^2 q + c12 p q^2 + c03 q^3 .
  !
  ! The functional derivative is implemented as the exact discrete adjoint of the
  ! discretized energy (transposes of the finite-difference first-derivative and
  ! radial-Laplacian stencils). The finite-difference derivative check can
  ! therefore be passed to high accuracy, and the ITE potential is strictly
  ! consistent with the evaluated energy.

  private
  public :: sl_is_active, sl_activate_manual, sl_deactivate, sl_load_coefficients, sl_report
  public :: sl_compute_pq, sl_F_and_derivs, sl_energy, sl_potential
  public :: sl_ncoeff, sl_basis_name, sl_get_params
  public :: apply_d1_transpose, apply_radlap_transpose

  integer, parameter :: sl_ncoeff = 8
  character(len=3), parameter :: sl_basis_name(sl_ncoeff) = &
       (/ 'c10', 'c20', 'c11', 'c02', 'c30', 'c21', 'c12', 'c03' /)

  logical  :: sl_active = .false.
  real(dp) :: sl_c(sl_ncoeff) = 0.0_dp
  real(dp) :: sl_pclip = 25.0_dp
  real(dp) :: sl_qclip = 25.0_dp
  real(dp) :: sl_rho_switch = 1.0e-4_dp
  real(dp) :: sl_rho_floor  = 1.0e-12_dp
  ! density-independent center regulator omega(r) = r^2/(r^2 + r_center^2),
  ! which removes the r->0 boundary artifact of the discrete Laplacian adjoint
  ! while keeping the functional exactly variational.
  real(dp) :: sl_r_center = 0.0_dp   ! 0 => omega(r)=1 (regulator off; opt-in via file)
  logical :: sl_pade = .false.        ! form: poly (default) or pade rational
  real(dp) :: sl_d(3) = 0.0_dp        ! Pade denominator: D = 1 + d10 p + d20 p^2 + d02 q^2

  real(dp), parameter :: a_pq = 4.0_dp * (3.0_dp*pi*pi)**(2.0_dp/3.0_dp)

contains

  logical function sl_is_active()
    sl_is_active = sl_active
  end function sl_is_active


  subroutine sl_deactivate()
    sl_active = .false.
  end subroutine sl_deactivate


  subroutine sl_activate_manual(c, pclip, qclip, rho_switch)
    real(dp), intent(in) :: c(sl_ncoeff), pclip, qclip, rho_switch
    sl_c(:) = c(:)
    sl_pclip = pclip
    sl_qclip = qclip
    sl_rho_switch = rho_switch
    sl_active = .true.
  end subroutine sl_activate_manual


  subroutine sl_get_params(c, pclip, qclip, rho_switch)
    real(dp), intent(out) :: c(sl_ncoeff), pclip, qclip, rho_switch
    c(:) = sl_c(:)
    pclip = sl_pclip
    qclip = sl_qclip
    rho_switch = sl_rho_switch
  end subroutine sl_get_params


  subroutine sl_load_coefficients(filename, ok)
    ! Read the coefficient file written by fit_semilocal.py:
    !   # comment lines
    !   pclip       <value>
    !   qclip       <value>
    !   rho_switch  <value>
    !   c10         <value>
    !   ...
    character(len=*), intent(in) :: filename
    logical, intent(out) :: ok

    integer :: unit, ios, k
    character(len=256) :: line
    character(len=32)  :: name
    real(dp) :: val
    logical :: found(sl_ncoeff)

    ok = .false.
    found(:) = .false.
    sl_c(:) = 0.0_dp

    unit = 71
    open(unit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(*,*) 'ERROR: cannot open semilocal coefficient file: ', trim(filename)
       return
    end if

    do
       read(unit,'(A)', iostat=ios) line
       if (ios /= 0) exit
       line = adjustl(line)
       if (len_trim(line) == 0) cycle
       if (line(1:1) == '#') cycle
       read(line, *, iostat=ios) name, val
       if (ios /= 0) then
          ! non-numeric value: only the 'form' key is allowed here
          block
            character(len=32) :: sval
            integer :: ios2
            read(line, *, iostat=ios2) name, sval
            if (ios2 == 0 .and. trim(name) == 'form') then
               sl_pade = (trim(sval) == 'pade')
            end if
          end block
          cycle
       end if

       if (trim(name) == 'form') then
          sl_pade = (val > 0.5_dp)
          cycle
       end if
       select case (trim(name))
       case ('pclip')
          sl_pclip = val
       case ('qclip')
          sl_qclip = val
       case ('rho_switch')
          sl_rho_switch = val
       case ('rho_floor_sl')
          sl_rho_floor = val
       case ('r_center')
          sl_r_center = val
       case ('d10', 'd1')
          sl_d(1) = val
       case ('d20', 'd2')
          sl_d(2) = val
       case ('d02', 'd3')
          sl_d(3) = val
       case default
          do k = 1, sl_ncoeff
             if (trim(name) == trim(sl_basis_name(k))) then
                sl_c(k) = val
                found(k) = .true.
             end if
          end do
       end select
    end do
    close(unit)

    if (.not. any(found)) then
       write(*,*) 'ERROR: no basis coefficients found in: ', trim(filename)
       return
    end if

    sl_active = .true.
    ok = .true.
    call sl_report()
  end subroutine sl_load_coefficients


  subroutine sl_report()
    integer :: k
    write(*,*) '-----------------------------------------------'
    write(*,*) 'Semilocal correction functional is ACTIVE'
    write(*,'(A,ES14.6,A,ES14.6)') '  pclip = ', sl_pclip, '  qclip = ', sl_qclip
    write(*,'(A,ES14.6,A,ES14.6)') '  rho_switch = ', sl_rho_switch, '  rho_floor_sl = ', sl_rho_floor
    write(*,'(A,ES14.6)') '  r_center = ', sl_r_center
    if (sl_pade) then
       write(*,'(A)') '  form = pade (rational, bounded)'
       write(*,'(A,3ES14.6)') '  d10 d20 d02 = ', sl_d(1), sl_d(2), sl_d(3)
    end if
    do k = 1, sl_ncoeff
       write(*,'(2X,A,A,ES18.10)') sl_basis_name(k), ' = ', sl_c(k)
    end do
    write(*,*) '-----------------------------------------------'
  end subroutine sl_report


  subroutine sl_compute_pq(nr, r, dr, rho, pvar, qvar, drho, laprho)
    ! Raw (unclipped) reduced invariants on the floored density.
    ! Used both by the fitting-dataset writer and internally.
    integer, intent(in) :: nr
    real(dp), intent(in) :: r(nr), dr, rho(nr)
    real(dp), intent(out) :: pvar(nr), qvar(nr), drho(nr), laprho(nr)

    real(dp), allocatable :: rho_f(:), d2rho(:)
    integer :: i

    allocate(rho_f(nr), d2rho(nr))
    do i = 1, nr
       rho_f(i) = max(rho(i), sl_rho_floor)
    end do

    call first_derivative_uniform(nr, dr, rho_f, drho)
    call second_derivative_uniform(nr, dr, rho_f, d2rho)
    do i = 1, nr
       laprho(i) = d2rho(i) + 2.0_dp*drho(i)/max(r(i), 1.0e-12_dp)
       pvar(i) = drho(i)**2 / (a_pq * rho_f(i)**(8.0_dp/3.0_dp))
       qvar(i) = laprho(i)  / (a_pq * rho_f(i)**(5.0_dp/3.0_dp))
    end do
    deallocate(rho_f, d2rho)
  end subroutine sl_compute_pq


  subroutine sl_F_and_derivs(p, q, F, Fp, Fq)
    ! F evaluated on the clipped invariants; Fp, Fq are the TOTAL derivatives
    ! dF/dp, dF/dq including the tanh clipping chain rule.
    real(dp), intent(in)  :: p, q
    real(dp), intent(out) :: F, Fp, Fq

    real(dp) :: tp, tq, pe, qe, dpe, dqe
    real(dp) :: Fpe, Fqe

    tp = tanh(p / sl_pclip)
    tq = tanh(q / sl_qclip)
    pe = sl_pclip * tp
    qe = sl_qclip * tq
    dpe = 1.0_dp - tp*tp
    dqe = 1.0_dp - tq*tq

    if (sl_pade) then
       ! F = N/D,  N = c10 pe + c20 pe^2 + c11 pe qe + c02 qe^2,
       !           D = 1 + d10 pe + d20 pe^2 + d02 qe^2  (d >= 0 => D >= 1).
       ! deg N = deg D  =>  F is globally bounded: no clip-saturation spikes
       ! in the Euler-Lagrange potential.
       block
         real(dp) :: Nv, Dv, Npe, Nqe, Dpe, Dqe
         Nv  = sl_c(1)*pe + sl_c(2)*pe*pe + sl_c(3)*pe*qe + sl_c(4)*qe*qe
         Dv  = 1.0_dp + sl_d(1)*pe + sl_d(2)*pe*pe + sl_d(3)*qe*qe
         Npe = sl_c(1) + 2.0_dp*sl_c(2)*pe + sl_c(3)*qe
         Nqe = sl_c(3)*pe + 2.0_dp*sl_c(4)*qe
         Dpe = sl_d(1) + 2.0_dp*sl_d(2)*pe
         Dqe = 2.0_dp*sl_d(3)*qe
         F   = Nv / Dv
         Fpe = (Npe*Dv - Nv*Dpe) / (Dv*Dv)
         Fqe = (Nqe*Dv - Nv*Dqe) / (Dv*Dv)
       end block
    else
       F = sl_c(1)*pe + sl_c(2)*pe*pe + sl_c(3)*pe*qe + sl_c(4)*qe*qe &
         + sl_c(5)*pe**3 + sl_c(6)*pe*pe*qe + sl_c(7)*pe*qe*qe + sl_c(8)*qe**3

       Fpe = sl_c(1) + 2.0_dp*sl_c(2)*pe + sl_c(3)*qe &
           + 3.0_dp*sl_c(5)*pe*pe + 2.0_dp*sl_c(6)*pe*qe + sl_c(7)*qe*qe
       Fqe = sl_c(3)*pe + 2.0_dp*sl_c(4)*qe &
           + sl_c(6)*pe*pe + 2.0_dp*sl_c(7)*pe*qe + 3.0_dp*sl_c(8)*qe*qe
    end if

    Fp = Fpe * dpe
    Fq = Fqe * dqe
  end subroutine sl_F_and_derivs


  real(dp) function sl_energy(nr, r, dr, w, rho)
    ! T_SL [MeV] for ONE species density.
    integer, intent(in) :: nr
    real(dp), intent(in) :: r(nr), dr, w(nr), rho(nr)

    real(dp), allocatable :: pvar(:), qvar(:), drho(:), laprho(:), edens(:)
    real(dp) :: rho_f, F, Fp, Fq, S
    integer :: i

    sl_energy = 0.0_dp
    if (.not. sl_active) return

    allocate(pvar(nr), qvar(nr), drho(nr), laprho(nr), edens(nr))
    call sl_compute_pq(nr, r, dr, rho, pvar, qvar, drho, laprho)

    do i = 1, nr
       rho_f = max(rho(i), sl_rho_floor)
       call sl_F_and_derivs(pvar(i), qvar(i), F, Fp, Fq)
       S = rho_f*rho_f / (rho_f*rho_f + sl_rho_switch*sl_rho_switch)
       S = S * r(i)*r(i) / (r(i)*r(i) + sl_r_center*sl_r_center)
       edens(i) = rho_f**(5.0_dp/3.0_dp) * S * F
    end do

    sl_energy = hbar2_over_2m * integrate_r3(nr, r, w, edens)
    deallocate(pvar, qvar, drho, laprho, edens)
  end function sl_energy


  subroutine sl_potential(nr, r, dr, w, rho, v_sl)
    ! Exact discrete gradient of the discretized T_SL:
    !   v_sl(j) = (dT_SL/drho_j) / V_j,   V_j = 4 pi r_j^2 w_j.
    ! Local part + D1^T and Lap^T adjoint contributions.
    integer, intent(in) :: nr
    real(dp), intent(in) :: r(nr), dr, w(nr), rho(nr)
    real(dp), intent(out) :: v_sl(nr)

    real(dp), allocatable :: pvar(:), qvar(:), drho(:), laprho(:)
    real(dp), allocatable :: alpha(:), beta(:), ta(:), tb(:), vol(:), mask(:)
    real(dp) :: rho_f, F, Fp, Fq, S, dS, denom
    integer :: i

    v_sl(:) = 0.0_dp
    if (.not. sl_active) return

    allocate(pvar(nr), qvar(nr), drho(nr), laprho(nr))
    allocate(alpha(nr), beta(nr), ta(nr), tb(nr), vol(nr), mask(nr))

    call sl_compute_pq(nr, r, dr, rho, pvar, qvar, drho, laprho)

    do i = 1, nr
       vol(i) = 4.0_dp*pi*r(i)**2 * w(i)
       if (rho(i) > sl_rho_floor) then
          mask(i) = 1.0_dp
       else
          mask(i) = 0.0_dp
       end if
       rho_f = max(rho(i), sl_rho_floor)
       call sl_F_and_derivs(pvar(i), qvar(i), F, Fp, Fq)
       denom = rho_f*rho_f + sl_rho_switch*sl_rho_switch
       S  = rho_f*rho_f / denom
       dS = 2.0_dp*rho_f*sl_rho_switch*sl_rho_switch / (denom*denom)
       ! density-independent center regulator (carried exactly through the adjoint)
       S  = S  * r(i)*r(i) / (r(i)*r(i) + sl_r_center*sl_r_center)
       dS = dS * r(i)*r(i) / (r(i)*r(i) + sl_r_center*sl_r_center)

       ! explicit-rho part of d e / d rho, e = rho^(5/3) S F(p,q)
       v_sl(i) = rho_f**(2.0_dp/3.0_dp) * S * &
            ( (5.0_dp/3.0_dp)*F - (8.0_dp/3.0_dp)*pvar(i)*Fp - (5.0_dp/3.0_dp)*qvar(i)*Fq ) &
            + rho_f**(5.0_dp/3.0_dp) * dS * F

       ! d e / d (D1 rho) and d e / d (Lap rho)
       alpha(i) = 2.0_dp * S * Fp * drho(i) / (a_pq * rho_f)
       beta(i)  = S * Fq / a_pq
    end do

    ! adjoint (transpose) contributions
    call apply_d1_transpose(nr, dr, vol(:)*alpha(:), ta)
    call apply_radlap_transpose(nr, r, dr, vol(:)*beta(:), tb)

    do i = 1, nr
       v_sl(i) = v_sl(i) + (ta(i) + tb(i)) / vol(i)
       v_sl(i) = hbar2_over_2m * v_sl(i) * mask(i)
    end do

    deallocate(pvar, qvar, drho, laprho, alpha, beta, ta, tb, vol, mask)
  end subroutine sl_potential


  subroutine apply_d1_transpose(n, dr, x, y)
    ! y = D1^T x, where D1 is the stencil of first_derivative_uniform:
    !   row 1:      (f2 - f1)/dr
    !   rows 2..n-1: (f_{i+1} - f_{i-1})/(2 dr)
    !   row n:      (f_n - f_{n-1})/dr
    integer, intent(in) :: n
    real(dp), intent(in) :: dr, x(n)
    real(dp), intent(out) :: y(n)
    integer :: i

    y(:) = 0.0_dp
    ! row 1
    y(1) = y(1) - x(1)/dr
    y(2) = y(2) + x(1)/dr
    ! interior rows
    do i = 2, n-1
       y(i-1) = y(i-1) - x(i)/(2.0_dp*dr)
       y(i+1) = y(i+1) + x(i)/(2.0_dp*dr)
    end do
    ! row n
    y(n-1) = y(n-1) - x(n)/dr
    y(n)   = y(n)   + x(n)/dr
  end subroutine apply_d1_transpose


  subroutine apply_radlap_transpose(n, r, dr, x, y)
    ! y = L^T x, where L is the stencil of radial_laplacian_uniform:
    !   L = D2 + diag(2/r) D1,
    ! with D2 from second_derivative_uniform:
    !   row 1:       (f3 - 2 f2 + f1)/dr^2
    !   rows 2..n-1: (f_{i+1} - 2 f_i + f_{i-1})/dr^2
    !   row n:       (f_n - 2 f_{n-1} + f_{n-2})/dr^2
    integer, intent(in) :: n
    real(dp), intent(in) :: r(n), dr, x(n)
    real(dp), intent(out) :: y(n)
    integer :: i
    real(dp) :: dr2, ci

    dr2 = dr*dr
    y(:) = 0.0_dp

    ! ---- D2^T part ----
    ! row 1
    y(1) = y(1) + x(1)/dr2
    y(2) = y(2) - 2.0_dp*x(1)/dr2
    y(3) = y(3) + x(1)/dr2
    do i = 2, n-1
       y(i-1) = y(i-1) + x(i)/dr2
       y(i)   = y(i)   - 2.0_dp*x(i)/dr2
       y(i+1) = y(i+1) + x(i)/dr2
    end do
    ! row n
    y(n-2) = y(n-2) + x(n)/dr2
    y(n-1) = y(n-1) - 2.0_dp*x(n)/dr2
    y(n)   = y(n)   + x(n)/dr2

    ! ---- (diag(2/r) D1)^T = D1^T diag(2/r) part ----
    ! row 1
    ci = 2.0_dp/max(r(1), 1.0e-12_dp) * x(1)
    y(1) = y(1) - ci/dr
    y(2) = y(2) + ci/dr
    do i = 2, n-1
       ci = 2.0_dp/max(r(i), 1.0e-12_dp) * x(i)
       y(i-1) = y(i-1) - ci/(2.0_dp*dr)
       y(i+1) = y(i+1) + ci/(2.0_dp*dr)
    end do
    ci = 2.0_dp/max(r(n), 1.0e-12_dp) * x(n)
    y(n-1) = y(n-1) - ci/dr
    y(n)   = y(n)   + ci/dr
  end subroutine apply_radlap_transpose

end module semilocal_correction_module

! The legacy module and namelist identifiers retain "teacher" for input-file
! compatibility; here the term denotes generation of reference fitting data.
module dd_teacher_module
  use precision_module
  use constants_module
  use radial_grid_module, only: compute_tf_vw_tau, integrate_r3
  use gauss_module
  use lindhard_kernel_module
  use nlf_postprocess_module, only: compute_constkf_convolution, compute_tau_nl_and_nlf, &
       write_kernel_file, write_nlf_file
  use semilocal_correction_module, only: sl_compute_pq
  implicit none

  ! kF-scale of the variational dd kernel: t = dd_kf_scale * 3 pi^2 (rho+rho')/2.
  ! The kernel normalization int f1 d^3R = 1 is kF-independent, so the TF limit
  ! is preserved exactly for ANY dd_kf_scale > 0.
  real(dp) :: dd_kf_scale = 1.0_dp
  ! pair-distance cache for the dd convolutions
  logical :: pd_ready = .false.
  integer :: pd_nr = -1, pd_nmu = -1
  real(dp) :: pd_rmax = -1.0_dp
  real(dp), allocatable :: pd_R(:,:,:), pd_wmu(:)
contains

  subroutine compute_ddkf_convolution(nr, r, w, rho, nmu, conv_dd)
    ! Density-dependent-kF reference convolution:
    !   I_dd(r) = int f1( kF_dd(r,r') |r-r'| ) rho^(5/6)(r') d^3r',
    !   kF_dd(r,r') = ( 3 pi^2 * ( rho(r)+rho(r') )/2 )^(1/3).
    ! Uses the scaling  f1(R; kF) = kF^3 * g(kF R)  with the dimensionless
    ! kernel table g(y) prepared by build_dimensionless_f1_table.
    integer, intent(in) :: nr, nmu
    real(dp), intent(in) :: r(nr), w(nr), rho(nr)
    real(dp), intent(out) :: conv_dd(nr)

    real(dp), allocatable :: mu(:), wmu(:), g56(:)
    integer :: i, j, k
    real(dp) :: ri, rj, Rij, angular_sum, kf3, kf, fval

    allocate(mu(nmu), wmu(nmu), g56(nr))
    call gauss_legendre(-1.0_dp, 1.0_dp, mu, wmu, nmu)

    do j = 1, nr
       g56(j) = max(rho(j), 0.0_dp)**(5.0_dp/6.0_dp)
    end do

    conv_dd(:) = 0.0_dp
    !$omp parallel do private(i, j, k, ri, rj, Rij, angular_sum, kf3, kf, fval) schedule(static)
    do i = 1, nr
       ri = r(i)
       do j = 1, nr
          rj = r(j)
          kf3 = 3.0_dp*pi*pi * 0.5_dp * (max(rho(i),0.0_dp) + max(rho(j),0.0_dp))
          if (kf3 <= 0.0_dp) cycle
          kf  = kf3**(1.0_dp/3.0_dp)
          angular_sum = 0.0_dp
          do k = 1, nmu
             Rij = sqrt(max(0.0_dp, ri*ri + rj*rj - 2.0_dp*ri*rj*mu(k)))
             fval = kf3 * interp_g_dimless(kf*Rij)
             angular_sum = angular_sum + wmu(k) * fval
          end do
          conv_dd(i) = conv_dd(i) + 2.0_dp*pi * rj*rj * w(j) * g56(j) * angular_sum
       end do
    end do
    !$omp end parallel do

    deallocate(mu, wmu, g56)
  end subroutine compute_ddkf_convolution


  subroutine compute_ddkf_conv_dconv(nr, r, w, rho, nmu, conv_dd, dconv_dd)
    ! Density-dependent-kF convolutions needed by the variational dd functional:
    !   conv_dd(i)  = sum_j V_j rho_j^(5/6) < f1(R; t_ij) >_mu ,
    !   dconv_dd(i) = sum_j V_j rho_j^(5/6) < d f1 / d t (R; t_ij) >_mu ,
    ! with t_ij = kF^3 = 3 pi^2 (rho_i + rho_j)/2 and the exact scaling
    !   f1(R; t) = t g(t^(1/3) R)  =>  d f1/d t = g(y) + (y/3) g'(y),  y = kF R.
    ! g and g' come from the same natural cubic spline of the dimensionless
    ! table, so the pair (energy, gradient) is exactly consistent.
    ! Exploits the i<->j symmetry of (t_ij, R_ij): each pair's angular sums are
    ! computed once and scattered to both rows (exact 2x saving); the pair
    ! distances R_ij(mu_k) are cached across calls since the grid is fixed.
    integer, intent(in) :: nr, nmu
    real(dp), intent(in) :: r(nr), w(nr), rho(nr)
    real(dp), intent(out) :: conv_dd(nr), dconv_dd(nr)

    real(dp), allocatable :: g56(:), vw56(:)
    integer :: i, j, k
    real(dp) :: kf3, kf, y, gv, gs, fsum, dsum

    call ensure_pairdist_cache(nr, r, nmu)

    allocate(g56(nr), vw56(nr))
    do j = 1, nr
       g56(j) = max(rho(j), 0.0_dp)**(5.0_dp/6.0_dp)
       vw56(j) = 2.0_dp*pi * r(j)*r(j) * w(j) * g56(j)
    end do

    conv_dd(:) = 0.0_dp
    dconv_dd(:) = 0.0_dp
    do i = 1, nr
       do j = i, nr
          kf3 = dd_kf_scale * 3.0_dp*pi*pi * 0.5_dp * (max(rho(i),0.0_dp) + max(rho(j),0.0_dp))
          if (kf3 <= 0.0_dp) cycle
          kf  = kf3**(1.0_dp/3.0_dp)
          fsum = 0.0_dp
          dsum = 0.0_dp
          do k = 1, nmu
             y = kf * pd_R(k, j, i)
             call interp_g_slope(y, gv, gs)
             fsum = fsum + pd_wmu(k) * kf3 * gv
             dsum = dsum + pd_wmu(k) * (gv + (y/3.0_dp)*gs)
          end do
          conv_dd(i)  = conv_dd(i)  + vw56(j) * fsum
          dconv_dd(i) = dconv_dd(i) + vw56(j) * dsum
          if (j /= i) then
             conv_dd(j)  = conv_dd(j)  + vw56(i) * fsum
             dconv_dd(j) = dconv_dd(j) + vw56(i) * dsum
          end if
       end do
    end do

    deallocate(g56, vw56)
  end subroutine compute_ddkf_conv_dconv


  subroutine ensure_pairdist_cache(nr, r, nmu)
    integer, intent(in) :: nr, nmu
    real(dp), intent(in) :: r(nr)
    integer :: i, j, k
    real(dp), allocatable :: mu(:)

    if (pd_ready .and. pd_nr == nr .and. pd_nmu == nmu) then
       if (abs(pd_rmax - r(nr)) < 1.0e-14_dp) return
    end if

    if (allocated(pd_R)) deallocate(pd_R)
    if (allocated(pd_wmu)) deallocate(pd_wmu)
    allocate(pd_R(nmu, nr, nr), pd_wmu(nmu), mu(nmu))
    call gauss_legendre(-1.0_dp, 1.0_dp, mu, pd_wmu, nmu)
    do i = 1, nr
       do j = i, nr
          do k = 1, nmu
             pd_R(k, j, i) = sqrt(max(0.0_dp, r(i)*r(i) + r(j)*r(j) - 2.0_dp*r(i)*r(j)*mu(k)))
          end do
       end do
    end do
    pd_ready = .true.
    pd_nr = nr
    pd_nmu = nmu
    pd_rmax = r(nr)
    deallocate(mu)
  end subroutine ensure_pairdist_cache


  subroutine run_dd_teacher(nr, r, dr, w, rho, kf_const, qmax, nq, nmu, nr_tab, &
       xqmax_dimless, ny_gtab, weight_alpha, teacher_file, dataset_file)
    ! Density-dependent-kF reference-data postprocessing on a given (KS) density:
    !  (1) constant-kF baseline tau_c and NLF_c,
    !  (2) density-dependent-kF reference tau_dd and NLF_dd,
    !  (3) fitting dataset:
    !      r, rho, drho, laprho, p, q, tau_c, tau_dd, y=(tau_dd-tau_c)/rho^(5/3), weight.
    integer, intent(in) :: nr, nq, nmu, nr_tab, ny_gtab
    real(dp), intent(in) :: r(nr), dr, w(nr), rho(nr)
    real(dp), intent(in) :: kf_const, qmax, xqmax_dimless, weight_alpha
    character(len=*), intent(in) :: teacher_file, dataset_file

    real(dp), allocatable :: rtab(:), ftab(:)
    real(dp), allocatable :: conv_c(:), conv_dd(:)
    real(dp), allocatable :: tau_tf(:), tau_vw(:), tau_c(:), nlf_c(:), tau_dd(:), nlf_dd(:)
    real(dp), allocatable :: pvar(:), qvar(:), drho(:), laprho(:)
    real(dp) :: rtab_max, kfmax, ymax, rho_f, yfit, wfit, dtau
    integer :: i, unit

    allocate(rtab(nr_tab), ftab(nr_tab))
    allocate(conv_c(nr), conv_dd(nr))
    allocate(tau_tf(nr), tau_vw(nr), tau_c(nr), nlf_c(nr), tau_dd(nr), nlf_dd(nr))
    allocate(pvar(nr), qvar(nr), drho(nr), laprho(nr))

    write(*,*) '-----------------------------------------------'
    write(*,*) 'Starting density-dependent-kF reference-data postprocessing'
    write(*,*) 'kf_const baseline [fm^-1] = ', kf_const
    write(*,*) 'weight_alpha              = ', weight_alpha
    write(*,*) '-----------------------------------------------'

    ! constant-kF baseline
    rtab_max = 2.0_dp * r(nr)
    call build_f1_table(kf_const, qmax, nq, nr_tab, rtab_max, rtab, ftab)
    call check_f1_norm(nr_tab, rtab, ftab)
    call compute_tf_vw_tau(nr, r, dr, rho, tau_tf, tau_vw)
    call compute_constkf_convolution(nr, r, w, rho, nmu, nr_tab, rtab, ftab, conv_c)
    call compute_tau_nl_and_nlf(nr, rho, conv_c, tau_tf, tau_vw, tau_c, nlf_c)

    ! dimensionless kernel table for the density-dependent kF
    kfmax = (3.0_dp*pi*pi*max(maxval(rho), 1.0e-8_dp))**(1.0_dp/3.0_dp)
    ymax  = 1.05_dp * kfmax * 2.0_dp * r(nr)
    call build_dimensionless_f1_table(xqmax_dimless, nq, ny_gtab, ymax)

    call compute_ddkf_convolution(nr, r, w, rho, nmu, conv_dd)
    call compute_tau_nl_and_nlf(nr, rho, conv_dd, tau_tf, tau_vw, tau_dd, nlf_dd)

    ! Semilocal invariants for the surrogate fit.
    call sl_compute_pq(nr, r, dr, rho, pvar, qvar, drho, laprho)

    ! Reference NLF comparison file.
    unit = 63
    open(unit=unit, file=trim(teacher_file), status='replace', action='write')
    write(unit,'(A)') '# density-dependent-kF reference-data postprocessing'
    write(unit,'(A)') '# columns: 1 r[fm]  2 rho_q  3 tau_TF  4 tau_vW  5 tau_NL_constkF  6 NLF_constkF  7 tau_NL_ddkF  8 NLF_ddkF  9 delta_tau=tau_dd-tau_c'
    do i = 1, nr
       write(unit,'(9ES24.15)') r(i), rho(i), tau_tf(i), tau_vw(i), tau_c(i), nlf_c(i), &
            tau_dd(i), nlf_dd(i), tau_dd(i)-tau_c(i)
    end do
    close(unit)

    ! Fitting dataset for the semilocal surrogate.
    unit = 64
    open(unit=unit, file=trim(dataset_file), status='replace', action='write')
    write(unit,'(A)') '# fitting dataset for the semilocal surrogate functional'
    write(unit,'(A,ES16.8)') '# weight_alpha = ', weight_alpha
    write(unit,'(A)') '# columns: 1 r[fm]  2 rho_q[fm^-3]  3 drho/dr  4 lap_rho  5 p=s^2  6 q  7 tau_constkF  8 tau_ddkF  9 y=(tau_dd-tau_c)/rho^(5/3)  10 weight=4*pi*r^2*dr*rho^alpha'
    do i = 1, nr
       rho_f = max(rho(i), 1.0e-12_dp)
       dtau  = tau_dd(i) - tau_c(i)
       yfit  = dtau / rho_f**(5.0_dp/3.0_dp)
       wfit  = 4.0_dp*pi*r(i)**2 * w(i) * rho_f**weight_alpha
       write(unit,'(10ES24.15)') r(i), rho(i), drho(i), laprho(i), pvar(i), qvar(i), &
            tau_c(i), tau_dd(i), yfit, wfit
    end do
    close(unit)

    write(*,*) 'Reference NLF output written to: ', trim(teacher_file)
    write(*,*) 'Fitting dataset written to     : ', trim(dataset_file)

    deallocate(rtab, ftab, conv_c, conv_dd)
    deallocate(tau_tf, tau_vw, tau_c, nlf_c, tau_dd, nlf_dd)
    deallocate(pvar, qvar, drho, laprho)
  end subroutine run_dd_teacher

  subroutine run_nlf_dd(nr, r, dr, w, rho, xqmax, nq, ny_base, nmu, nlf_file)
    ! Evaluate the density-dependent kernel, including any loaded delta-g
    ! correction, through the shared dimensionless g table. The reported NLF
    ! therefore uses the same functional representation as a dd-mode SCF run.
    integer, intent(in) :: nr, nq, ny_base, nmu
    real(dp), intent(in) :: r(nr), dr, w(nr), rho(nr), xqmax
    character(len=*), intent(in) :: nlf_file

    real(dp), allocatable :: conv(:), tau_tf(:), tau_vw(:), tau_nl(:), nlf(:)
    integer :: ny_gtab
    real(dp) :: ymax

    allocate(conv(nr), tau_tf(nr), tau_vw(nr), tau_nl(nr), nlf(nr))
    ! Use the same dimensionless-table rule as the SCF/derivative path:
    ! adaptive ymax/ny scaled to the radial box, with the
    ! namelist-controlled xqmax, nq, and ny base, so the NLF postprocessor and
    ! the SCF/derivative paths share one explicitly recorded discretization.
    ymax = 1.05_dp * (3.0_dp*pi*pi*0.5_dp)**(1.0_dp/3.0_dp) * 2.0_dp * r(nr)
    ny_gtab = max(ny_base, int(ymax*125.0_dp) + 1)
    call build_dimensionless_f1_table(xqmax, nq, ny_gtab, ymax)

    write(*,'(A,F6.1,A,I0,A,I0,A,F7.2)') ' Starting dd-kernel NLF postprocessing'// &
         ' (delta-g aware): xqmax=', xqmax, ' nq=', nq, &
         ' ny=', ny_gtab, ' ymax=', ymax
    call compute_tf_vw_tau(nr, r, dr, rho, tau_tf, tau_vw)
    ! Use the same convolution routine as the SCF/ITE density evolution
    ! (compute_ddkf_conv_dconv), so the NLF diagnostic is consistent with the
    ! functional that generated the density. We keep only conv.
    block
      real(dp), allocatable :: dconv_tmp(:)
      allocate(dconv_tmp(nr))
      call compute_ddkf_conv_dconv(nr, r, w, rho, nmu, conv, dconv_tmp)
      deallocate(dconv_tmp)
    end block

    call compute_tau_nl_and_nlf(nr, rho, conv, tau_tf, tau_vw, tau_nl, nlf)
    call write_nlf_file(nlf_file, nr, r, rho, conv, tau_tf, tau_vw, tau_nl, nlf)
    write(*,*) 'dd NLF output written to: ', trim(nlf_file)
    deallocate(conv, tau_tf, tau_vw, tau_nl, nlf)
  end subroutine run_nlf_dd

end module dd_teacher_module


module nl_ofdft_derivative_module
  use precision_module
  use constants_module
  use radial_grid_module, only: radial_laplacian_uniform, first_derivative_uniform, &
       compute_tf_vw_tau, integrate_r3
  use skyrme_skp_module, only: compute_skp_potential
  use lindhard_kernel_module, only: build_f1_table, check_f1_norm, build_dimensionless_f1_table
  use nlf_postprocess_module, only: compute_constkf_convolution, compute_tau_nl_and_nlf, &
       write_kernel_file
  use dd_teacher_module, only: compute_ddkf_conv_dconv, dd_kf_scale, run_nlf_dd
  use semilocal_correction_module, only: sl_is_active, sl_potential, sl_energy
  implicit none

  ! --- kernel mode switches (set from the main program) ---
  logical :: dd_kernel_mode = .false.
  real(dp) :: dd_xqmax = 80.0_dp
  integer :: dd_nq = 6000
  integer :: dd_ny_base = 8001
contains

  subroutine run_ofdft_derivative_constkf(nr, r, dr, w, rho_q, &
       kf_const, qmax, nq, nmu, nr_tab, derivative_file, derivative_energy_file, kernel_file)
    integer, intent(in) :: nr, nq, nmu, nr_tab
    real(dp), intent(in) :: r(nr), dr, w(nr), rho_q(nr)
    real(dp), intent(in) :: kf_const, qmax
    character(len=*), intent(in) :: derivative_file, derivative_energy_file, kernel_file

    real(dp), allocatable :: rtab(:), ftab(:), conv(:), dconv(:)
    real(dp), allocatable :: tau_tf(:), tau_vw(:), tau_nl(:), nlf_nl(:)
    real(dp), allocatable :: phi(:), lap_phi(:), rho_tot(:), v_skp(:)
    real(dp), allocatable :: v_vw(:), v_kernel(:), v_corr(:), v_sl(:)
    real(dp), allocatable :: dTdrho(:), dEdrho(:), v_eff_ite(:), hphi(:), res_phi(:)
    real(dp) :: rtab_max, mu_est
    real(dp) :: T_vw, T_nl, T_sl, E_skp, Etot_like, norm_q, residual_norm
    real(dp) :: E_grad

    allocate(rtab(nr_tab), ftab(nr_tab), conv(nr), dconv(nr))
    allocate(tau_tf(nr), tau_vw(nr), tau_nl(nr), nlf_nl(nr))
    allocate(phi(nr), lap_phi(nr), rho_tot(nr), v_skp(nr))
    allocate(v_vw(nr), v_kernel(nr), v_corr(nr), v_sl(nr))
    allocate(dTdrho(nr), dEdrho(nr), v_eff_ite(nr), hphi(nr), res_phi(nr))

    rtab_max = 2.0_dp * r(nr)

    write(*,*) '-----------------------------------------------'
    write(*,*) 'Starting NL-KEDF functional-derivative evaluation'
    write(*,*) 'reference kf_const [fm^-1] = ', kf_const
    write(*,*) 'qmax [fm^-1]     = ', qmax
    write(*,*) 'nq, nmu, nr_tab  = ', nq, nmu, nr_tab
    write(*,*) '-----------------------------------------------'

    call build_f1_table(kf_const, qmax, nq, nr_tab, rtab_max, rtab, ftab)
    call check_f1_norm(nr_tab, rtab, ftab)
    call write_kernel_file(kernel_file, nr_tab, rtab, ftab)

    call compute_kernel_conv(nr, r, w, rho_q, nmu, nr_tab, rtab, ftab, conv, dconv)
    call compute_tf_vw_tau(nr, r, dr, rho_q, tau_tf, tau_vw)
    call compute_tau_nl_and_nlf(nr, rho_q, conv, tau_tf, tau_vw, tau_nl, nlf_nl)

    rho_tot(:) = 2.0_dp * rho_q(:)
    call compute_skp_potential(nr, r, dr, rho_tot, v_skp)

    call compute_constkf_derivative_arrays(nr, r, dr, w, rho_q, conv, dconv, v_skp, &
         phi, lap_phi, v_vw, v_kernel, v_corr, v_sl, dTdrho, dEdrho, v_eff_ite, hphi, &
         mu_est, res_phi, residual_norm)

    T_vw = hbar2_over_2m * integrate_r3(nr, r, w, tau_vw)
    T_nl = hbar2_over_2m * integrate_r3(nr, r, w, tau_nl)
    T_sl = sl_energy(nr, r, dr, w, rho_q)
    E_skp = compute_skp_energy(nr, r, dr, w, rho_tot)
    ! Spatial-gradient regularization energy, per species; consistent
    ! with the EL potential in compute_constkf_derivative_arrays and with
    ! evaluate_ite_energies. E_grad_species = s2 (hbar^2/2m) int rho^(1/3)|grad rho|^2.
    E_grad = 0.0_dp
    if (abs(s2_grad) > 0.0_dp) then
       block
         real(dp), allocatable :: rp(:), eg(:)
         integer :: ig
         allocate(rp(nr), eg(nr))
         call first_derivative_uniform(nr, dr, rho_q, rp)
         do ig = 1, nr
            eg(ig) = s2_grad * max(rho_q(ig),rho_floor)**(1.0_dp/3.0_dp) * rp(ig)**2
         end do
         E_grad = hbar2_over_2m * integrate_r3(nr, r, w, eg)
         deallocate(rp, eg)
       end block
    end if
    Etot_like = 2.0_dp*(T_nl + T_sl + E_grad) + E_skp
    norm_q = integrate_r3(nr, r, w, rho_q)

    call write_derivative_radial_file(derivative_file, nr, r, rho_q, rho_tot, conv, phi, lap_phi, &
         tau_tf, tau_vw, tau_nl, nlf_nl, v_vw, v_kernel, v_corr, v_sl, dTdrho, v_skp, &
         dEdrho, v_eff_ite, hphi, res_phi, mu_est)

    call write_derivative_energy_file(derivative_energy_file, norm_q, T_vw, T_nl, T_sl, E_skp, &
         Etot_like, mu_est, residual_norm, E_grad)

    write(*,*) 'Derivative radial output written to: ', trim(derivative_file)
    write(*,*) 'Derivative summary written to      : ', trim(derivative_energy_file)

    deallocate(rtab, ftab, conv, tau_tf, tau_vw, tau_nl, nlf_nl)
    deallocate(phi, lap_phi, rho_tot, v_skp, v_vw, v_kernel, v_corr, v_sl)
    deallocate(dTdrho, dEdrho, v_eff_ite, hphi, res_phi)
  end subroutine run_ofdft_derivative_constkf


  subroutine compute_kernel_conv(nr, r, w, rho, nmu, nr_tab, rtab, ftab, conv, dconv)
    ! Kernel convolution dispatcher.  kernel_mode 'const': constant-kF table
    ! convolution, dconv = 0.  kernel_mode 'dd': exactly-variational
    ! density-dependent-kF convolutions (see compute_ddkf_conv_dconv).
    integer, intent(in) :: nr, nmu, nr_tab
    real(dp), intent(in) :: r(nr), w(nr), rho(nr)
    real(dp), intent(in) :: rtab(nr_tab), ftab(nr_tab)
    real(dp), intent(out) :: conv(nr), dconv(nr)
    real(dp) :: ymax_need
    integer :: ny_need

    if (dd_kernel_mode) then
       ! Generous cap: species densities never approach 0.5 fm^-3.
       ymax_need = 1.05_dp * (3.0_dp*pi*pi*0.5_dp)**(1.0_dp/3.0_dp) * 2.0_dp * r(nr)
       ny_need = max(dd_ny_base, int(ymax_need*125.0_dp) + 1)
       call build_dimensionless_f1_table(dd_xqmax, dd_nq, ny_need, ymax_need)
       call compute_ddkf_conv_dconv(nr, r, w, rho, nmu, conv, dconv)
    else
       call compute_constkf_convolution(nr, r, w, rho, nmu, nr_tab, rtab, ftab, conv)
       dconv(:) = 0.0_dp
    end if
  end subroutine compute_kernel_conv


  subroutine compute_constkf_derivative_arrays(nr, r, dr, w, rho_q, conv, dconv, v_skp, &
       phi, lap_phi, v_vw, v_kernel, v_corr, v_sl, dTdrho, dEdrho, v_eff_ite, hphi, &
       mu_est, res_phi, residual_norm)
    integer, intent(in) :: nr
    real(dp), intent(in) :: r(nr), dr, w(nr), rho_q(nr), conv(nr), dconv(nr), v_skp(nr)
    real(dp), intent(out) :: phi(nr), lap_phi(nr)
    real(dp), intent(out) :: v_vw(nr), v_kernel(nr), v_corr(nr), v_sl(nr)
    real(dp), intent(out) :: dTdrho(nr), dEdrho(nr), v_eff_ite(nr), hphi(nr)
    real(dp), intent(out) :: mu_est, res_phi(nr), residual_norm

    integer :: i
    real(dp) :: rho_safe, denom, norm_phi, mu_num, res2
    real(dp), allocatable :: rho_p(:), rho_pp(:), v_grad2(:)

    do i = 1, nr
       phi(i) = sqrt(max(rho_q(i), rho_floor))
    end do

    call radial_laplacian_uniform(nr, r, dr, phi, lap_phi)

    ! Spatial-gradient regularization potential with low-density suppression.
    ! Exact discrete adjoint of E3 = sum_i W_i rho_i^(1/3) (D rho)_i^2:
    !   dE3/drho_j = W_j (1/3) rho_j^(-2/3) (Drho_j)^2 + (D^T [ W 2 rho^(1/3) Drho ])_j,
    ! divided by the volume weight to give a per-d^3r potential. Finite in the tail.
    allocate(rho_p(nr), rho_pp(nr), v_grad2(nr))
    v_grad2(:) = 0.0_dp
    if (abs(s2_grad) > 0.0_dp) then
       block
         real(dp), allocatable :: Drho(:), gvec(:), DTg(:), Wvol(:)
         integer :: j
         allocate(Drho(nr), gvec(nr), DTg(nr), Wvol(nr))
         do j = 1, nr
            Wvol(j) = 4.0_dp*pi*r(j)*r(j)*w(j)
         end do
         call first_derivative_uniform(nr, dr, rho_q, Drho)
         do j = 1, nr
            gvec(j) = 2.0_dp * Wvol(j) * max(rho_q(j),rho_floor)**(1.0_dp/3.0_dp) * Drho(j)
         end do
         DTg(:) = 0.0_dp
         do j = 2, nr-1
            DTg(j+1) = DTg(j+1) + gvec(j)/(2.0_dp*dr)
            DTg(j-1) = DTg(j-1) - gvec(j)/(2.0_dp*dr)
         end do
         DTg(2)    = DTg(2)    + gvec(1)/dr
         DTg(1)    = DTg(1)    - gvec(1)/dr
         DTg(nr)   = DTg(nr)   + gvec(nr)/dr
         DTg(nr-1) = DTg(nr-1) - gvec(nr)/dr
         do j = 1, nr
            v_grad2(j) = s2_grad * hbar2_over_2m * &
                 ( Wvol(j)*(1.0_dp/3.0_dp)*max(rho_q(j),rho_floor)**(-2.0_dp/3.0_dp)*Drho(j)**2 &
                   + DTg(j) ) / max(Wvol(j), 1.0e-30_dp)
         end do
         deallocate(Drho, gvec, DTg, Wvol)
       end block
    end if

    ! Fitted semilocal correction potential (zero when the fit is inactive).
    if (sl_is_active()) then
       call sl_potential(nr, r, dr, w, rho_q, v_sl)
    else
       v_sl(:) = 0.0_dp
    end if

    do i = 1, nr
       rho_safe = max(rho_q(i), rho_floor)

       ! d T_vW / d rho = -C laplacian(sqrt(rho)) / sqrt(rho).
       v_vw(i) = -hbar2_over_2m * lap_phi(i) / phi(i)

       ! Functional derivative of the two-point nonlocal term.  In dd mode the
       ! second piece is the kernel's own density dependence:
       !   d t_ij / d rho_x = (3 pi^2 / 2)(delta_ix + delta_jx)  (symmetric),
       !   => v_dd = C c_NL [ (5/3) rho^(-1/6) conv + 3 pi^2 rho^(5/6) dconv ].
       v_kernel(i) = hbar2_over_2m * c_NL * ( (5.0_dp/3.0_dp) * &
            rho_safe**(-1.0_dp/6.0_dp) * conv(i) &
            + dd_kf_scale * 3.0_dp*pi*pi * max(rho_q(i),0.0_dp)**(5.0_dp/6.0_dp) * dconv(i) )

       ! Functional derivative of the local -B int rho^(5/3) correction, plus
       ! the spatial-gradient regularization potential.
       v_corr(i) = -hbar2_over_2m * (5.0_dp/3.0_dp) * c_CORR * &
            rho_safe**(2.0_dp/3.0_dp) + v_grad2(i)

       dTdrho(i) = v_vw(i) + v_kernel(i) + v_corr(i) + v_sl(i)
       dEdrho(i) = dTdrho(i) + v_skp(i)

       ! Local potential part for a Schrodinger-like ITE implementation.
       ! The vW derivative is represented by the kinetic operator -C laplacian.
       v_eff_ite(i) = v_kernel(i) + v_corr(i) + v_sl(i) + v_skp(i)
       hphi(i) = -hbar2_over_2m * lap_phi(i) + v_eff_ite(i) * phi(i)
    end do
    deallocate(rho_p, rho_pp, v_grad2)

    norm_phi = 0.0_dp
    mu_num = 0.0_dp
    do i = 1, nr
       norm_phi = norm_phi + 4.0_dp*pi*r(i)**2 * phi(i)**2 * w(i)
       mu_num = mu_num + 4.0_dp*pi*r(i)**2 * phi(i) * hphi(i) * w(i)
    end do

    if (norm_phi > 0.0_dp) then
       mu_est = mu_num / norm_phi
    else
       mu_est = 0.0_dp
    end if

    res2 = 0.0_dp
    do i = 1, nr
       res_phi(i) = hphi(i) - mu_est * phi(i)
       res2 = res2 + 4.0_dp*pi*r(i)**2 * res_phi(i)**2 * w(i)
    end do
    residual_norm = sqrt(max(res2, 0.0_dp))
  end subroutine compute_constkf_derivative_arrays


  real(dp) function compute_skp_energy(nr, r, dr, w, rho_tot)
    integer, intent(in) :: nr
    real(dp), intent(in) :: r(nr), dr, w(nr), rho_tot(nr)

    real(dp), allocatable :: drho(:), edens(:)
    integer :: i
    real(dp) :: cgrad

    allocate(drho(nr), edens(nr))

    call first_derivative_uniform(nr, dr, rho_tot, drho)
    cgrad = 9.0_dp*t1_skp - 5.0_dp*t2_skp - 4.0_dp*t2_skp*x2_skp

    do i = 1, nr
       edens(i) = (3.0_dp/8.0_dp)*t0_skp*rho_tot(i)**2 &
            + (1.0_dp/16.0_dp)*t3_skp*rho_tot(i)**(2.0_dp + gamma_skp) &
            + (1.0_dp/64.0_dp)*cgrad*drho(i)**2
    end do

    compute_skp_energy = integrate_r3(nr, r, w, edens)
    deallocate(drho, edens)
  end function compute_skp_energy


  subroutine write_derivative_radial_file(filename, nr, r, rho_q, rho_tot, conv, phi, lap_phi, &
       tau_tf, tau_vw, tau_nl, nlf_nl, v_vw, v_kernel, v_corr, v_sl, dTdrho, v_skp, &
       dEdrho, v_eff_ite, hphi, res_phi, mu_est)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nr
    real(dp), intent(in) :: r(nr), rho_q(nr), rho_tot(nr), conv(nr), phi(nr), lap_phi(nr)
    real(dp), intent(in) :: tau_tf(nr), tau_vw(nr), tau_nl(nr), nlf_nl(nr)
    real(dp), intent(in) :: v_vw(nr), v_kernel(nr), v_corr(nr), v_sl(nr), dTdrho(nr), v_skp(nr)
    real(dp), intent(in) :: dEdrho(nr), v_eff_ite(nr), hphi(nr), res_phi(nr), mu_est

    integer :: i, unit

    unit = 61
    open(unit=unit, file=trim(filename), status='replace', action='write')
    write(unit,'(A)') '# Nonlocal-kernel NL-OFDFT functional-derivative radial output'
    write(unit,'(A)') '# rho_q is one-species density; rho_tot = 2*rho_q for N=Z SkP potential.'
    write(unit,'(A,ES24.15)') '# mu_est[MeV] = ', mu_est
    write(unit,'(A)') '# columns:'
    write(unit,'(A)') '#  1 r[fm]'
    write(unit,'(A)') '#  2 rho_q[fm^-3]'
    write(unit,'(A)') '#  3 rho_tot[fm^-3]'
    write(unit,'(A)') '#  4 conv = int f1(kF|r-rp|) rho_q(rp)^(5/6) d^3rp'
    write(unit,'(A)') '#  5 phi=sqrt(rho_q)'
    write(unit,'(A)') '#  6 lap_phi[fm^(-7/2)]'
    write(unit,'(A)') '#  7 tau_TF'
    write(unit,'(A)') '#  8 tau_vW'
    write(unit,'(A)') '#  9 tau_NL'
    write(unit,'(A)') '# 10 NLF_NL'
    write(unit,'(A)') '# 11 v_vW=dT_vW/drho[MeV] = -C lap_phi/phi'
    write(unit,'(A)') '# 12 v_kernel=dT_nonlocal_kernel/drho[MeV]'
    write(unit,'(A)') '# 13 v_corr=dT_local_correction/drho[MeV]'
    write(unit,'(A)') '# 14 dT_NL/drho[MeV] = v_vW+v_kernel+v_corr'
    write(unit,'(A)') '# 15 v_Skyrme=dE_int/drho_tot[MeV]'
    write(unit,'(A)') '# 16 dE_total/drho_q[MeV] = dT_NL/drho_q+v_Skyrme'
    write(unit,'(A)') '# 17 v_eff_ITE[MeV] = v_kernel+v_corr+v_Skyrme, excluding vW kinetic operator'
    write(unit,'(A)') '# 18 H_OF phi'
    write(unit,'(A)') '# 19 residual_phi = H_OF phi - mu_est phi'
    write(unit,'(A)') '# 20 v_SL = d T_SL^fit / drho[MeV] (zero if the semilocal fit is inactive)'

    do i = 1, nr
       write(unit,'(20ES24.15)') r(i), rho_q(i), rho_tot(i), conv(i), phi(i), lap_phi(i), &
            tau_tf(i), tau_vw(i), tau_nl(i), nlf_nl(i), v_vw(i), v_kernel(i), v_corr(i), &
            dTdrho(i), v_skp(i), dEdrho(i), v_eff_ite(i), hphi(i), res_phi(i), v_sl(i)
    end do
    close(unit)
  end subroutine write_derivative_radial_file


  subroutine write_derivative_energy_file(filename, norm_q, T_vw, T_nl, T_sl, E_skp, Etot_like, &
       mu_est, residual_norm, E_grad)
    character(len=*), intent(in) :: filename
    real(dp), intent(in) :: norm_q, T_vw, T_nl, T_sl, E_skp, Etot_like, mu_est, residual_norm
    real(dp), intent(in) :: E_grad

    integer :: unit
    unit = 62
    open(unit=unit, file=trim(filename), status='replace', action='write')
    write(unit,'(A)') '# Nonlocal-kernel NL-OFDFT functional-derivative summary'
    write(unit,'(A,ES24.15)') 'N_species              ', norm_q
    write(unit,'(A,ES24.15)') 'T_vW_species_MeV       ', T_vw
    ! T_KEDFnl includes vW + nonlocal kernel + subtractive local rho^(5/3) correction
    ! (it is the full tau_nl integral, not the pure two-point term).
    write(unit,'(A,ES24.15)') 'T_KEDFnl_species_MeV    ', T_nl
    write(unit,'(A,ES24.15)') 'E_grad_species_MeV      ', E_grad
    write(unit,'(A,ES24.15)') 'T_SL_fit_species_MeV    ', T_sl
    write(unit,'(A,ES24.15)') 'E_SkP_total_MeV         ', E_skp
    write(unit,'(A,ES24.15)') 'E_total_NeqZ_MeV        ', Etot_like
    write(unit,'(A,ES24.15)') 'mu_est_MeV             ', mu_est
    write(unit,'(A,ES24.15)') 'residual_phi_L2        ', residual_norm
    close(unit)
  end subroutine write_derivative_energy_file

end module nl_ofdft_derivative_module

module derivative_checker_module
  use precision_module
  use constants_module
  use radial_grid_module, only: first_derivative_uniform, integrate_r3
  use gauss_module
  use lindhard_kernel_module
  use nlf_postprocess_module, only: compute_constkf_convolution
  use nl_ofdft_derivative_module, only: compute_kernel_conv
  use dd_teacher_module, only: dd_kf_scale
  use semilocal_correction_module, only: sl_is_active, sl_energy, sl_potential, &
       apply_d1_transpose
  implicit none

  ! Finite-difference derivative checker.
  !
  ! The checked energy is the one the code actually varies:
  !   E_sym[rho_q] = 2 T_model[rho_q] + E_int[2 rho_q],
  !   T_model = T_vW + T_kernel + T_corr (+ T_SL^fit if active),
  ! where T_kernel is selected by kernel_mode.
  ! The analytic directional derivative is
  !   D_ana = sum_j V_j * (dE_sym/drho_q)_j * eta_j,   V_j = 4 pi r_j^2 w_j,
  ! computed from exact discrete gradients (adjoint stencils) of the discretized
  ! energies, and is compared with
  !   D_num(eps) = ( E[rho + eps eta] - E[rho - eps eta] ) / (2 eps)
  ! for eps = 1e-2 ... 1e-7.  The perturbation is multiplicative,
  !   eta(r) = rho(r) * exp( -((r-r0)/sigma)^2 ),
  ! so that rho +- eps eta stays positive and no density floor is activated.

contains

  subroutine run_derivative_checker(nr, r, dr, w, rho, kf_const, qmax, nq, nmu, nr_tab, &
       eta_r0, eta_sigma, out_file)
    integer, intent(in) :: nr, nq, nmu, nr_tab
    real(dp), intent(in) :: r(nr), dr, w(nr), rho(nr), kf_const, qmax
    real(dp), intent(in) :: eta_r0, eta_sigma
    character(len=*), intent(in) :: out_file

    integer, parameter :: ncomp = 6   ! vW, kernel, corr, SL, SkP, total
    integer, parameter :: neps = 11
    character(len=8), parameter :: comp_name(ncomp) = &
         (/ 'vW      ', 'kernel  ', 'corr    ', 'SL_fit  ', 'SkP     ', 'TOTAL   ' /)

    real(dp), allocatable :: rtab(:), ftab(:), eta(:), vol(:)
    real(dp), allocatable :: grad(:,:), rho_p(:), rho_m(:)
    real(dp) :: rtab_max, eps, lg
    real(dp) :: E_p(ncomp), E_m(ncomp), Dana(ncomp), Dnum(ncomp), relerr(ncomp)
    real(dp) :: maxrel(ncomp)
    integer :: i, ie, ic, unit

    allocate(rtab(nr_tab), ftab(nr_tab), eta(nr), vol(nr))
    allocate(grad(nr, ncomp), rho_p(nr), rho_m(nr))

    write(*,*) '-----------------------------------------------'
    write(*,*) 'Starting finite-difference derivative checker'
    write(*,'(A,F8.3,A,F8.3)') ' eta bump: r0 = ', eta_r0, '  sigma = ', eta_sigma
    if (sl_is_active()) then
       write(*,*) ' semilocal correction: INCLUDED in the checked energy'
    else
       write(*,*) ' semilocal correction: not active'
    end if
    write(*,*) '-----------------------------------------------'

    rtab_max = 2.0_dp * r(nr)
    call build_f1_table(kf_const, qmax, nq, nr_tab, rtab_max, rtab, ftab)

    do i = 1, nr
       vol(i) = 4.0_dp*pi*r(i)**2 * w(i)
       eta(i) = rho(i) * exp(-((r(i)-eta_r0)/eta_sigma)**2)
    end do

    ! ---- analytic discrete gradients dE/drho_q at rho ----
    call assemble_discrete_gradients(nr, r, dr, w, rho, nmu, nr_tab, rtab, ftab, ncomp, grad)

    do ic = 1, ncomp
       Dana(ic) = 0.0_dp
       do i = 1, nr
          Dana(ic) = Dana(ic) + vol(i) * grad(i, ic) * eta(i)
       end do
    end do

    unit = 65
    open(unit=unit, file=trim(out_file), status='replace', action='write')
    write(unit,'(A)') '# finite-difference derivative checker for E_sym[rho_q] = 2 T_model + E_int[2 rho_q]'
    write(unit,'(A)') '# eta(r) = rho(r) * exp(-((r-r0)/sigma)^2)'
    write(unit,'(A,F10.4,A,F10.4)') '# r0 = ', eta_r0, '  sigma = ', eta_sigma
    write(unit,'(A,6ES20.10)') '# D_ana per component (vW kernel corr SL SkP TOTAL): ', (Dana(ic), ic=1,ncomp)
    write(unit,'(A)') '# columns: 1 eps  2 Dnum_TOTAL  3 Dana_TOTAL  4 relerr_TOTAL  5 relerr_vW  6 relerr_kernel  7 relerr_corr  8 relerr_SL  9 relerr_SkP'

    maxrel(:) = 0.0_dp
    do ie = 1, neps
       lg = -2.0_dp - 5.0_dp*real(ie-1,dp)/real(neps-1,dp)   ! log10(eps): -2 ... -7
       eps = 10.0_dp**lg

       rho_p(:) = rho(:) + eps*eta(:)
       rho_m(:) = rho(:) - eps*eta(:)

       call evaluate_component_energies(nr, r, dr, w, rho_p, nmu, nr_tab, rtab, ftab, ncomp, E_p)
       call evaluate_component_energies(nr, r, dr, w, rho_m, nmu, nr_tab, rtab, ftab, ncomp, E_m)

       do ic = 1, ncomp
          Dnum(ic) = (E_p(ic) - E_m(ic)) / (2.0_dp*eps)
          if (abs(Dana(ic)) > 1.0e-14_dp) then
             relerr(ic) = abs(Dnum(ic) - Dana(ic)) / abs(Dana(ic))
          else
             relerr(ic) = abs(Dnum(ic) - Dana(ic))
          end if
          maxrel(ic) = max(maxrel(ic), relerr(ic))
       end do

       write(unit,'(9ES18.8)') eps, Dnum(ncomp), Dana(ncomp), relerr(ncomp), &
            relerr(1), relerr(2), relerr(3), relerr(4), relerr(5)
       write(*,'(A,ES10.2,A,ES12.4,A,ES12.4)') '  eps=', eps, '  relerr(TOTAL)=', relerr(ncomp), &
            '  Dnum(TOTAL)=', Dnum(ncomp)
    end do

    write(unit,'(A)') '#'
    write(unit,'(A)') '# component-resolved minimum relative error over the eps sweep is reported on screen.'
    close(unit)

    write(*,*) 'Derivative checker output written to: ', trim(out_file)
    do ic = 1, ncomp
       write(*,'(2A,ES12.4)') '  max relerr over sweep, ', comp_name(ic), maxrel(ic)
    end do

    deallocate(rtab, ftab, eta, vol, grad, rho_p, rho_m)
  end subroutine run_derivative_checker


  subroutine evaluate_component_energies(nr, r, dr, w, rho, nmu, nr_tab, rtab, ftab, ncomp, E)
    ! E(1) = 2 T_vW, E(2) = 2 T_kernel, E(3) = 2 T_corr, E(4) = 2 T_SL,
    ! E(5) = E_SkP[2 rho], E(6) = total E_sym.
    integer, intent(in) :: nr, nmu, nr_tab, ncomp
    real(dp), intent(in) :: r(nr), dr, w(nr), rho(nr)
    real(dp), intent(in) :: rtab(nr_tab), ftab(nr_tab)
    real(dp), intent(out) :: E(ncomp)

    real(dp), allocatable :: sq(:), dsq(:), conv(:), dconv(:), edens(:), rho_tot(:), drt(:)
    real(dp) :: cgrad
    integer :: i

    allocate(sq(nr), dsq(nr), conv(nr), dconv(nr), edens(nr), rho_tot(nr), drt(nr))

    ! vW
    do i = 1, nr
       sq(i) = sqrt(max(rho(i), 0.0_dp))
    end do
    call first_derivative_uniform(nr, dr, sq, dsq)
    do i = 1, nr
       edens(i) = dsq(i)**2
    end do
    E(1) = 2.0_dp * hbar2_over_2m * integrate_r3(nr, r, w, edens)

    ! kernel (const or dd depending on kernel_mode; the energy density has the
    ! same form  c_NL rho^(5/6) conv  in both cases)
    call compute_kernel_conv(nr, r, w, rho, nmu, nr_tab, rtab, ftab, conv, dconv)
    do i = 1, nr
       edens(i) = c_NL * max(rho(i),0.0_dp)**(5.0_dp/6.0_dp) * conv(i)
    end do
    E(2) = 2.0_dp * hbar2_over_2m * integrate_r3(nr, r, w, edens)

    ! local correction
    do i = 1, nr
       edens(i) = -c_CORR * max(rho(i),0.0_dp)**(5.0_dp/3.0_dp)
    end do
    E(3) = 2.0_dp * hbar2_over_2m * integrate_r3(nr, r, w, edens)

    ! Spatial-gradient regularization with low-density suppression:
    !   E3 = s2_grad * (hbar^2/2m) * int rho^(1/3) |grad rho|^2 d^3r.
    ! The rho^(1/3) prefactor removes the rho^(-4/3) singularity of the plain
    ! Weizsacker-like form, so the potential stays finite in the exponential
    ! tail; it vanishes for uniform matter, preserving the Thomas-Fermi limit.
    if (abs(s2_grad) > 0.0_dp) then
       block
         real(dp), allocatable :: rp(:)
         allocate(rp(nr))
         call first_derivative_uniform(nr, dr, rho, rp)
         do i = 1, nr
            edens(i) = s2_grad * max(rho(i),rho_floor)**(1.0_dp/3.0_dp) * rp(i)**2
         end do
         E(3) = E(3) + 2.0_dp * hbar2_over_2m * integrate_r3(nr, r, w, edens)
         deallocate(rp)
       end block
    end if

    ! fitted semilocal correction
    E(4) = 2.0_dp * sl_energy(nr, r, dr, w, rho)

    ! SkP interaction energy at rho_tot = 2 rho
    rho_tot(:) = 2.0_dp*rho(:)
    call first_derivative_uniform(nr, dr, rho_tot, drt)
    cgrad = 9.0_dp*t1_skp - 5.0_dp*t2_skp - 4.0_dp*t2_skp*x2_skp
    do i = 1, nr
       edens(i) = (3.0_dp/8.0_dp)*t0_skp*rho_tot(i)**2 &
            + (1.0_dp/16.0_dp)*t3_skp*max(rho_tot(i),0.0_dp)**(2.0_dp + gamma_skp) &
            + (1.0_dp/64.0_dp)*cgrad*drt(i)**2
    end do
    E(5) = integrate_r3(nr, r, w, edens)

    E(6) = E(1) + E(2) + E(3) + E(4) + E(5)

    deallocate(sq, dsq, conv, dconv, edens, rho_tot, drt)
  end subroutine evaluate_component_energies


  subroutine assemble_discrete_gradients(nr, r, dr, w, rho, nmu, nr_tab, rtab, ftab, ncomp, grad)
    ! Exact discrete gradients (dE/drho_q)_j / V_j of the component energies above.
    integer, intent(in) :: nr, nmu, nr_tab, ncomp
    real(dp), intent(in) :: r(nr), dr, w(nr), rho(nr)
    real(dp), intent(in) :: rtab(nr_tab), ftab(nr_tab)
    real(dp), intent(out) :: grad(nr, ncomp)

    real(dp), allocatable :: sq(:), dsq(:), conv(:), dconv(:), vol(:), tmp(:), ty(:)
    real(dp), allocatable :: rho_tot(:), drt(:)
    real(dp) :: cgrad, rh
    integer :: i

    allocate(sq(nr), dsq(nr), conv(nr), dconv(nr), vol(nr), tmp(nr), ty(nr))
    allocate(rho_tot(nr), drt(nr))

    do i = 1, nr
       vol(i) = 4.0_dp*pi*r(i)**2 * w(i)
    end do

    ! ---- vW: T = C sum V (D1 sqrt(rho))^2 ;
    !      dT/drho_j = C [D1^T (2 V D1 sqrt(rho))]_j / (2 sqrt(rho_j))
    do i = 1, nr
       sq(i) = sqrt(max(rho(i), rho_floor))
    end do
    call first_derivative_uniform(nr, dr, sq, dsq)
    tmp(:) = 2.0_dp * vol(:) * dsq(:)
    call apply_d1_transpose(nr, dr, tmp, ty)
    do i = 1, nr
       grad(i,1) = 2.0_dp * hbar2_over_2m * ty(i) / (2.0_dp*sq(i)) / vol(i)
    end do

    ! ---- kernel: const mode  dT/drho = C c_NL (5/3) rho^(-1/6) conv  (dconv = 0);
    !      dd mode adds the kernel's own density dependence
    !      + C c_NL 3 pi^2 rho^(5/6) dconv  (exact discrete gradient, see text).
    call compute_kernel_conv(nr, r, w, rho, nmu, nr_tab, rtab, ftab, conv, dconv)
    do i = 1, nr
       rh = max(rho(i), rho_floor)
       grad(i,2) = 2.0_dp * hbar2_over_2m * c_NL * ( (5.0_dp/3.0_dp) * rh**(-1.0_dp/6.0_dp) * conv(i) &
            + dd_kf_scale * 3.0_dp*pi*pi * max(rho(i),0.0_dp)**(5.0_dp/6.0_dp) * dconv(i) )
    end do

    ! ---- local correction plus spatial-gradient regularization
    do i = 1, nr
       rh = max(rho(i), 0.0_dp)
       grad(i,3) = -2.0_dp * hbar2_over_2m * (5.0_dp/3.0_dp) * c_CORR * rh**(2.0_dp/3.0_dp)
    end do
    if (abs(s2_grad) > 0.0_dp) then
       ! E3 = 2 * s2_grad * (hbar^2/2m) * sum V rho^(1/3) (D1 rho)^2.
       ! Discrete adjoint:
       !   dE3/drho_j = 2 s2 C [ V_j (1/3) rho_j^-2/3 (D1rho_j)^2
       !                         + (D1^T [ V 2 rho^1/3 D1rho ])_j ] / V_j.
       block
         real(dp), allocatable :: d1rho(:), gg(:), gt(:)
         integer :: j
         allocate(d1rho(nr), gg(nr), gt(nr))
         call first_derivative_uniform(nr, dr, rho, d1rho)
         do j = 1, nr
            gg(j) = 2.0_dp * vol(j) * max(rho(j),1.0e-20_dp)**(1.0_dp/3.0_dp) * d1rho(j)
         end do
         call apply_d1_transpose(nr, dr, gg, gt)
         do j = 1, nr
            grad(j,3) = grad(j,3) + 2.0_dp * s2_grad * hbar2_over_2m * &
                 ( vol(j)*(1.0_dp/3.0_dp)*max(rho(j),1.0e-20_dp)**(-2.0_dp/3.0_dp)*d1rho(j)**2 &
                   + gt(j) ) / max(vol(j),1.0e-30_dp)
         end do
         deallocate(d1rho, gg, gt)
       end block
    end if

    ! ---- fitted semilocal correction (exact discrete adjoint inside sl_potential)
    call sl_potential(nr, r, dr, w, rho, tmp)
    grad(:,4) = 2.0_dp * tmp(:)

    ! ---- SkP: E = sum V [ 3/8 t0 rt^2 + 1/16 t3 rt^(2+g) + 1/64 cgrad (D1 rt)^2 ],
    !      dE/drho_q = 2 * dE/drho_tot
    rho_tot(:) = 2.0_dp*rho(:)
    call first_derivative_uniform(nr, dr, rho_tot, drt)
    cgrad = 9.0_dp*t1_skp - 5.0_dp*t2_skp - 4.0_dp*t2_skp*x2_skp
    tmp(:) = (1.0_dp/64.0_dp)*cgrad * 2.0_dp * vol(:) * drt(:)
    call apply_d1_transpose(nr, dr, tmp, ty)
    do i = 1, nr
       rh = max(rho_tot(i), 0.0_dp)
       grad(i,5) = 2.0_dp * ( (3.0_dp/4.0_dp)*t0_skp*rho_tot(i) &
            + ((2.0_dp + gamma_skp)/16.0_dp)*t3_skp*rh**(1.0_dp + gamma_skp) &
            + ty(i)/vol(i) )
    end do

    grad(:,6) = grad(:,1) + grad(:,2) + grad(:,3) + grad(:,4) + grad(:,5)

    deallocate(sq, dsq, conv, dconv, vol, tmp, ty, rho_tot, drt)
  end subroutine assemble_discrete_gradients

end module derivative_checker_module




module direct_minimization_module
  ! Direct energy minimization (DFTpy / PROFESS style) as an alternative to the
  ! split-step imaginary-time evolution.  Minimizes E_sym[rho] over the
  ! norm-invariant parametrization
  !     rho_i = xNpart * phi_i^2 / <phi^2>,   <phi^2> = sum_j V_j phi_j^2,
  ! with the exact discrete gradient
  !     dE/dphi_i = 2 phi_i V_i c ( v_i - mubar ),
  !     c = xNpart/<phi^2>,  v_i = (dE/drho)_i / V_i,  mubar = (sum_j V_j v_j rho_j)/xNpart,
  ! built from the FD-verified gradient assembly of derivative_checker_module.
  ! Solver: Polak-Ribiere+ nonlinear conjugate gradient with Armijo backtracking
  ! line search (monotone by construction, NO time-step / CFL restriction).
  use precision_module
  use constants_module
  use lindhard_kernel_module, only: build_f1_table, check_f1_norm
  use derivative_checker_module, only: evaluate_component_energies, assemble_discrete_gradients
  implicit none
contains

  subroutine run_direct_minimization(nr, r, dr, w, rho_init, xNpart, kf_const, qmax, nq, nmu, nr_tab, &
       max_iter, gtol, history_file, density_file, summary_file)
    integer, intent(in) :: nr, nq, nmu, nr_tab, max_iter
    real(dp), intent(in) :: r(nr), dr, w(nr), rho_init(nr), xNpart, kf_const, qmax, gtol
    character(len=*), intent(in) :: history_file, density_file, summary_file

    integer, parameter :: ncomp = 6
    real(dp), allocatable :: rtab(:), ftab(:), vol(:)
    real(dp), allocatable :: phi(:), rho(:), g(:), d(:), g_old(:), phi_t(:), rho_t(:)
    real(dp), allocatable :: grad(:,:)
    real(dp) :: E(ncomp), E_t(ncomp)
    real(dp) :: rtab_max, gg, gg_old, gy, beta, alpha, alpha0, dginit, E0, mubar, c, s2
    real(dp) :: gnorm, dE_last
    integer :: i, it, ls, unit, n_ls, restarts
    logical :: ok

    allocate(rtab(nr_tab), ftab(nr_tab), vol(nr))
    allocate(phi(nr), rho(nr), g(nr), d(nr), g_old(nr), phi_t(nr), rho_t(nr), grad(nr, ncomp))

    rtab_max = 2.0_dp * r(nr)
    call build_f1_table(kf_const, qmax, nq, nr_tab, rtab_max, rtab, ftab)
    call check_f1_norm(nr_tab, rtab, ftab)

    do i = 1, nr
       vol(i) = 4.0_dp*pi*r(i)*r(i)*w(i)
       phi(i) = sqrt(max(rho_init(i), 0.0_dp))
    end do
    call renorm_rho(nr, vol, xNpart, phi, rho, c, s2)

    open(newunit=unit, file=history_file, status='replace', action='write')
    write(unit,'(A)') '# iter  E_sym  gnorm_phi  alpha  n_linesearch  E_vW2 E_ker2 E_corr2 E_SL2 E_SkP'

    call evaluate_component_energies(nr, r, dr, w, rho, nmu, nr_tab, rtab, ftab, ncomp, E)
    call gradient_phi(nr, r, dr, w, rho, vol, phi, c, xNpart, nmu, nr_tab, rtab, ftab, ncomp, grad, g, mubar)
    gg = dot_v(nr, g, g)
    gnorm = sqrt(gg)
    d = -g
    alpha0 = 1.0e-4_dp
    dE_last = huge(1.0_dp)
    restarts = 0

    write(*,'(A)') ' Direct minimization (PR+ CG, Armijo backtracking):'
    write(*,'(A,ES14.6,A,ES12.4)') '   E0 = ', E(6), '   |g|0 = ', gnorm

    do it = 1, max_iter
       dginit = dot_v(nr, g, d)
       if (dginit >= 0.0_dp) then
          d = -g
          dginit = -gg
          restarts = restarts + 1
       end if
       E0 = E(6)
       alpha = alpha0
       ok = .false.
       n_ls = 0
       do ls = 1, 40
          n_ls = n_ls + 1
          phi_t = phi + alpha*d
          call renorm_rho(nr, vol, xNpart, phi_t, rho_t, c, s2)
          call evaluate_component_energies(nr, r, dr, w, rho_t, nmu, nr_tab, rtab, ftab, ncomp, E_t)
          if (E_t(6) <= E0 + 1.0e-4_dp*alpha*dginit) then
             ok = .true.
             exit
          end if
          alpha = 0.5_dp*alpha
       end do
       if (.not. ok) then
          write(*,'(A,I8)') '   line search failed (converged or pathological) at iter ', it
          exit
       end if

       phi = phi_t
       rho = rho_t
       dE_last = E0 - E_t(6)
       E = E_t

       g_old = g
       gg_old = gg
       call gradient_phi(nr, r, dr, w, rho, vol, phi, c, xNpart, nmu, nr_tab, rtab, ftab, ncomp, grad, g, mubar)
       gg = dot_v(nr, g, g)
       gnorm = sqrt(gg)

       gy = gg - dot_v(nr, g, g_old)
       beta = max(0.0_dp, gy/max(gg_old, 1.0e-300_dp))
       if (mod(it, 4*nr) == 0) beta = 0.0_dp
       do i = 1, nr
          d(i) = -g(i) + beta*d(i)
       end do

       if (n_ls == 1) then
          alpha0 = 2.0_dp*alpha
       else
          alpha0 = alpha
       end if

       write(unit,'(I8,3ES20.10,I4,5ES16.8)') it, E(6), gnorm, alpha, n_ls, E(1), E(2), E(3), E(4), E(5)
       if (gnorm < gtol .and. dE_last < 1.0e-10_dp) exit
    end do

    close(unit)

    open(newunit=unit, file=density_file, status='replace', action='write')
    write(unit,'(A)') '# r[fm] rho_q_MIN[fm^-3]'
    do i = 1, nr
       write(unit,'(2ES24.15)') r(i), rho(i)
    end do
    close(unit)

    open(newunit=unit, file=summary_file, status='replace', action='write')
    write(unit,'(A,ES24.15)') 'E_total_NeqZ_MeV        ', E(6)
    write(unit,'(A,ES24.15)') 'mu_est_MeV              ', mubar
    write(unit,'(A,ES24.15)') 'gnorm_phi               ', gnorm
    write(unit,'(A,ES24.15)') 'T_vW2_MeV               ', E(1)
    write(unit,'(A,ES24.15)') 'T_kernel2_MeV           ', E(2)
    write(unit,'(A,ES24.15)') 'T_corr2_MeV             ', E(3)
    write(unit,'(A,ES24.15)') 'T_SL2_MeV               ', E(4)
    write(unit,'(A,ES24.15)') 'E_SkP_MeV               ', E(5)
    write(unit,'(A,I12)')     'cg_restarts             ', restarts
    close(unit)

    write(*,'(A,ES16.8,A,ES12.4)') ' Direct minimization done: E = ', E(6), '  |g| = ', gnorm
    write(*,'(A,A)') ' Minimized density written to: ', trim(density_file)

    deallocate(rtab, ftab, vol, phi, rho, g, d, g_old, phi_t, rho_t, grad)
  end subroutine run_direct_minimization


  subroutine renorm_rho(nr, vol, xNpart, phi, rho, c, s2)
    integer, intent(in) :: nr
    real(dp), intent(in) :: vol(nr), xNpart, phi(nr)
    real(dp), intent(out) :: rho(nr), c, s2
    integer :: i
    s2 = 0.0_dp
    do i = 1, nr
       s2 = s2 + vol(i)*phi(i)*phi(i)
    end do
    c = xNpart / max(s2, 1.0e-300_dp)
    do i = 1, nr
       rho(i) = c * phi(i)*phi(i)
    end do
  end subroutine renorm_rho


  subroutine gradient_phi(nr, r, dr, w, rho, vol, phi, c, xNpart, nmu, nr_tab, rtab, ftab, ncomp, grad, g, mubar)
    integer, intent(in) :: nr, nmu, nr_tab, ncomp
    real(dp), intent(in) :: r(nr), dr, w(nr), rho(nr), vol(nr), phi(nr), c, xNpart
    real(dp), intent(in) :: rtab(nr_tab), ftab(nr_tab)
    real(dp), intent(inout) :: grad(nr, ncomp)
    real(dp), intent(out) :: g(nr), mubar
    integer :: i
    call assemble_discrete_gradients(nr, r, dr, w, rho, nmu, nr_tab, rtab, ftab, ncomp, grad)
    mubar = 0.0_dp
    do i = 1, nr
       mubar = mubar + vol(i)*grad(i, ncomp)*rho(i)
    end do
    mubar = mubar / xNpart
    do i = 1, nr
       g(i) = 2.0_dp*phi(i)*vol(i)*c*(grad(i, ncomp) - mubar)
    end do
  end subroutine gradient_phi


  real(dp) function dot_v(n, a, b)
    integer, intent(in) :: n
    real(dp), intent(in) :: a(n), b(n)
    integer :: i
    dot_v = 0.0_dp
    do i = 1, n
       dot_v = dot_v + a(i)*b(i)
    end do
  end function dot_v

end module direct_minimization_module




module radial_initial_density_module
  use precision_module
  use constants_module, only: pi, rho_floor
  use radial_grid_module, only: integrate_r3
  implicit none
contains

  subroutine init_species_fermi_density(nr, r, w, A_mass, Nq, r0_nuc, a_surf, rho_q)
    integer, intent(in) :: nr
    real(dp), intent(in) :: r(nr), w(nr)
    real(dp), intent(in) :: A_mass, Nq, r0_nuc, a_surf
    real(dp), intent(out) :: rho_q(nr)

    integer :: i
    real(dp) :: R0, a0, x, shape, norm

    if (r0_nuc > 0.0_dp) then
       R0 = r0_nuc * A_mass**(1.0_dp/3.0_dp)
    else
       R0 = 1.20_dp * A_mass**(1.0_dp/3.0_dp)
    end if

    if (a_surf > 0.0_dp) then
       a0 = a_surf
    else
       a0 = 0.55_dp
    end if

    do i = 1, nr
       x = (r(i) - R0) / a0
       if (x > 80.0_dp) then
          shape = 0.0_dp
       else if (x < -80.0_dp) then
          shape = 1.0_dp
       else
          shape = 1.0_dp / (1.0_dp + exp(x))
       end if
       rho_q(i) = shape
    end do

    call renormalize_species_density(nr, r, w, Nq, rho_q)
  end subroutine init_species_fermi_density


  subroutine renormalize_species_density(nr, r, w, Nq, rho_q)
    integer, intent(in) :: nr
    real(dp), intent(in) :: r(nr), w(nr), Nq
    real(dp), intent(inout) :: rho_q(nr)
    real(dp) :: norm
    integer :: i

    do i = 1, nr
       rho_q(i) = max(rho_q(i), rho_floor)
    end do

    norm = integrate_r3(nr, r, w, rho_q)
    if (norm <= 0.0_dp) stop 'ERROR: density normalization failed in renormalize_species_density.'
    rho_q(:) = rho_q(:) * Nq / norm
  end subroutine renormalize_species_density


  subroutine density_to_phi(nr, rho_q, phi)
    integer, intent(in) :: nr
    real(dp), intent(in) :: rho_q(nr)
    real(dp), intent(out) :: phi(nr)
    integer :: i

    do i = 1, nr
       phi(i) = sqrt(max(rho_q(i), rho_floor))
    end do
  end subroutine density_to_phi


  subroutine phi_to_density(nr, phi, rho_q)
    integer, intent(in) :: nr
    real(dp), intent(in) :: phi(nr)
    real(dp), intent(out) :: rho_q(nr)
    integer :: i

    do i = 1, nr
       rho_q(i) = max(phi(i)*phi(i), rho_floor)
    end do
  end subroutine phi_to_density

end module radial_initial_density_module

module radial_ite_split_module
  use precision_module
  use constants_module
  use radial_grid_module, only: integrate_r3, radial_laplacian_uniform, first_derivative_uniform, &
       compute_tf_vw_tau
  use radial_initial_density_module, only: renormalize_species_density, density_to_phi, phi_to_density
  use lapack_tridiag_module, only: diagonalize_tridiagonal
  use skyrme_skp_module, only: compute_skp_potential
  use lindhard_kernel_module, only: build_f1_table, check_f1_norm
  use nlf_postprocess_module, only: compute_constkf_convolution, write_kernel_file
  use nl_ofdft_derivative_module, only: compute_kernel_conv, compute_constkf_derivative_arrays, &
       run_ofdft_derivative_constkf, compute_skp_energy
  use semilocal_correction_module, only: sl_is_active, sl_energy
  implicit none

  private
  public :: run_ite_constkf_radial
  public :: compute_ite_vloc_constkf, evaluate_ite_energies

  logical :: kinetic_ready = .false.
  integer :: kinetic_nr_saved = -1
  real(dp) :: kinetic_dr_saved = -1.0_dp
  real(dp), allocatable :: kinetic_eval(:), kinetic_evec(:,:)

contains

  subroutine run_ite_constkf_radial(nr, r, dr, w, rho_init, npart_species, method, dt, max_iter, tol, &
       fpi_max_iter, fpi_tol, kf_const, qmax, nq, nmu, nr_tab, monitor_every, &
       density_file, history_file, derivative_file, summary_file, kernel_file, &
       dt_adapt, dt_max, dt_grow, dt_shrink)
    integer, intent(in) :: nr, max_iter, fpi_max_iter, nq, nmu, nr_tab, monitor_every
    real(dp), intent(in) :: r(nr), dr, w(nr), rho_init(nr), npart_species
    real(dp), intent(in) :: dt, tol, fpi_tol, kf_const, qmax
    character(len=*), intent(in) :: method
    character(len=*), intent(in) :: density_file, history_file, derivative_file, summary_file, kernel_file
    logical, intent(in) :: dt_adapt
    real(dp), intent(in) :: dt_max, dt_grow, dt_shrink

    real(dp), allocatable :: rtab(:), ftab(:)
    real(dp), allocatable :: phi(:), phi_new(:), rho(:), rho_new(:)
    real(dp) :: rtab_max, delta, norm_now, maxrho
    real(dp) :: T_vw_m, T_ker_m, T_corr_m, T_sl_m, E_skp_m, E_sym_m, mu_m, resid_m
    integer :: iter, unit
    real(dp) :: dt_cur, E_prev, E_now
    real(dp), allocatable :: phi_prev(:)
    logical :: have_prev
    integer :: n_reject, n_accept, last_accepted_iter
    integer :: monitored_flag, accepted_flag, rollback_flag

    allocate(rtab(nr_tab), ftab(nr_tab))
    allocate(phi(nr), phi_new(nr), rho(nr), rho_new(nr))

    rho(:) = rho_init(:)
    call renormalize_species_density(nr, r, w, npart_species, rho)
    call density_to_phi(nr, rho, phi)
    call normalize_phi(nr, r, w, npart_species, phi)

    rtab_max = 2.0_dp * r(nr)
    write(*,*) '-----------------------------------------------'
    write(*,*) 'Starting radial NL-KEDF imaginary-time evolution'
    write(*,*) 'method          = ', trim(method)
    write(*,*) 'dt [MeV^-1]     = ', dt
    write(*,*) 'N_species       = ', npart_species
    write(*,*) 'reference kf_const [fm^-1]= ', kf_const
    write(*,*) '-----------------------------------------------'

    call build_f1_table(kf_const, qmax, nq, nr_tab, rtab_max, rtab, ftab)
    call check_f1_norm(nr_tab, rtab, ftab)
    call write_kernel_file(kernel_file, nr_tab, rtab, ftab)
    call init_radial_kinetic_propagator(nr, dr)

    T_vw_m = 0.0_dp; T_ker_m = 0.0_dp; T_corr_m = 0.0_dp; T_sl_m = 0.0_dp
    E_skp_m = 0.0_dp; E_sym_m = 0.0_dp; mu_m = 0.0_dp; resid_m = 0.0_dp

    unit = 81
    open(unit=unit, file=trim(history_file), status='replace', action='write')
    write(unit,'(A)') '# ITE history for the variational nonlocal-kernel NL-OFDFT (+ fitted semilocal correction if active)'
    write(unit,'(A)') '# energies are recomputed every monitor_every iterations; otherwise the last monitored value is repeated.'
    write(unit,'(A)') '# columns: 1 iter  2 delta_phi_L2  3 N_species  4 max_rho  5 E_sym[MeV]  6 T_vW_sp  7 T_kernel_sp  8 T_corr_sp(incl E_grad)  9 T_SL_sp  10 E_SkP  11 mu_est[MeV]  12 residual_phi_L2  13 dt  14 monitored  15 accepted  16 rollback'

    allocate(phi_prev(nr))
    dt_cur = dt
    have_prev = .false.
    n_reject = 0; n_accept = 0; last_accepted_iter = 0
    E_prev = 0.0_dp

    do iter = 1, max_iter
       call ite_step_constkf(nr, r, dr, w, phi, npart_species, method, dt_cur, fpi_max_iter, fpi_tol, &
            nmu, nr_tab, rtab, ftab, phi_new)

       delta = phi_l2_difference(nr, r, w, phi, phi_new)
       call phi_to_density(nr, phi_new, rho_new)
       norm_now = integrate_r3(nr, r, w, rho_new)
       maxrho = maxval(rho_new)

       monitored_flag = 0; accepted_flag = 1; rollback_flag = 0
       if (iter == 1 .or. mod(iter, max(1, monitor_every)) == 0) then
          monitored_flag = 1
          call evaluate_ite_energies(nr, r, dr, w, rho_new, nmu, nr_tab, rtab, ftab, &
               T_vw_m, T_ker_m, T_corr_m, T_sl_m, E_skp_m, E_sym_m, mu_m, resid_m)

          if (dt_adapt) then
             ! Adaptive imaginary-time step: on a monotone monitor interval grow
             ! dt (up to dt_max); on an energy increase, shrink dt and roll the
             ! state back to the previous monitored point, KEEPING that point as
             ! the reference so the retry must still improve on it. This enforces
             ! monotonic energy at the monitor points while annealing dt upward.
             E_now = E_sym_m
             if (have_prev) then
                if (E_now > E_prev + 1.0e-12_dp) then
                   ! reject: shrink dt, restore last good state, retry.
                   ! phi_prev/E_prev remain the reference (do NOT clear have_prev).
                   dt_cur = max(dt_cur*dt_shrink, 1.0e-8_dp)
                   phi(:) = phi_prev(:)
                   call phi_to_density(nr, phi, rho)
                   n_reject = n_reject + 1
                   rollback_flag = 1; accepted_flag = 0
                   write(unit,'(I8,11ES24.15,ES24.15,3I6)') iter, delta, norm_now, maxrho, &
                        E_sym_m, T_vw_m, T_ker_m, T_corr_m, T_sl_m, E_skp_m, mu_m, resid_m, &
                        dt_cur, monitored_flag, accepted_flag, rollback_flag
                   cycle
                else
                   dt_cur = min(dt_cur*dt_grow, dt_max)
                end if
             end if
             phi_prev(:) = phi_new(:)
             E_prev = E_now
             have_prev = .true.
             n_accept = n_accept + 1
             last_accepted_iter = iter
          end if
       end if

       write(unit,'(I8,11ES24.15,ES24.15,3I6)') iter, delta, norm_now, maxrho, E_sym_m, &
            T_vw_m, T_ker_m, T_corr_m, T_sl_m, E_skp_m, mu_m, resid_m, &
            dt_cur, monitored_flag, accepted_flag, rollback_flag

       phi(:) = phi_new(:)
       rho(:) = rho_new(:)

       if (mod(iter, max(1, max_iter/10)) == 0) then
          write(*,'(A,I8,A,ES12.4,A,ES12.4,A,F14.6,A,ES12.4)') 'ITE iter ', iter, ' delta=', delta, &
               ' N=', norm_now, ' E_sym=', E_sym_m, ' res=', resid_m
       end if
       if (delta < tol) exit
    end do
    ! Record attempted versus accepted iteration accounting.
    write(unit,'(A)') '# adaptive accounting:'
    write(unit,'(A,I0)') '# attempted_iterations ', min(iter, max_iter)
    write(unit,'(A,I0)') '# accepted_monitor_points ', n_accept
    write(unit,'(A,I0)') '# rollback_events ', n_reject
    write(unit,'(A,I0)') '# last_accepted_iteration ', last_accepted_iter
    close(unit)
    deallocate(phi_prev)

    call phi_to_density(nr, phi, rho)
    call renormalize_species_density(nr, r, w, npart_species, rho)
    call write_density_output(density_file, nr, r, rho)

    ! Final derivative/residual output at the ITE density.
    call run_ofdft_derivative_constkf(nr, r, dr, w, rho, kf_const, qmax, nq, nmu, nr_tab, &
         derivative_file, summary_file, kernel_file)

    write(*,*) 'ITE density output written to   : ', trim(density_file)
    write(*,*) 'ITE history output written to   : ', trim(history_file)
    write(*,*) 'ITE derivative output written to: ', trim(derivative_file)

    deallocate(rtab, ftab, phi, phi_new, rho, rho_new)
  end subroutine run_ite_constkf_radial


  subroutine ite_step_constkf(nr, r, dr, w, phi_in, npart_species, method, dt, fpi_max_iter, fpi_tol, &
       nmu, nr_tab, rtab, ftab, phi_out)
    integer, intent(in) :: nr, fpi_max_iter, nmu, nr_tab
    real(dp), intent(in) :: r(nr), dr, w(nr), phi_in(nr), npart_species, dt, fpi_tol
    real(dp), intent(in) :: rtab(nr_tab), ftab(nr_tab)
    character(len=*), intent(in) :: method
    real(dp), intent(out) :: phi_out(nr)

    real(dp), allocatable :: rho(:), vloc(:), tmp(:), tmp2(:), vgrad(:), vtilde(:)
    character(len=32) :: m

    allocate(rho(nr), vloc(nr), tmp(nr), tmp2(nr), vgrad(nr), vtilde(nr))
    m = adjustl(method)

    call phi_to_density(nr, phi_in, rho)
    call compute_ite_vloc_constkf(nr, r, dr, w, rho, nmu, nr_tab, rtab, ftab, vloc)

    select case (trim(m))
    case ('ITE_1A')
       call potential_step(nr, dt, 1.0_dp, vloc, phi_in, tmp)
       call kinetic_step_l0(nr, r, dt, 1.0_dp, tmp, phi_out)

    case ('ITE_1B0')
       call kinetic_step_l0(nr, r, dt, 1.0_dp, phi_in, tmp)
       call potential_step(nr, dt, 1.0_dp, vloc, tmp, phi_out)

    case ('ITE_1BF')
       call kinetic_step_l0(nr, r, dt, 1.0_dp, phi_in, tmp)
       call potential_fpi(nr, r, dr, w, tmp, phi_in, npart_species, dt, 1.0_dp, &
            fpi_max_iter, fpi_tol, nmu, nr_tab, rtab, ftab, phi_out)

    case ('ITE_2A0')
       call potential_step(nr, dt, 0.5_dp, vloc, phi_in, tmp)
       call kinetic_step_l0(nr, r, dt, 1.0_dp, tmp, tmp2)
       call potential_step(nr, dt, 0.5_dp, vloc, tmp2, phi_out)

    case ('ITE_SI2')
       ! Semi-implicit second-order (Strang) step: kinetic half-step (exact
       ! exponential of the frozen -D2 operator, unconditionally stable),
       ! then a full implicit potential step (backward-Euler on the local
       ! potential with a few inner fixed-point sweeps so vloc[rho] stays
       ! consistent), then a kinetic half-step. Second-order accurate and
       ! implicit in BOTH operators, so it keeps the large-dt stability of
       ! ITE_BE while removing its first-order splitting error.
       call kinetic_step_l0(nr, r, dt, 0.5_dp, phi_in, tmp)
       call potential_implicit_step(nr, r, dr, w, dt, tmp, nmu, nr_tab, rtab, ftab, tmp2)
       call kinetic_step_l0(nr, r, dt, 0.5_dp, tmp2, phi_out)

    case ('ITE_BE')
       ! Fully-implicit backward-Euler step on u=r*phi:
       !   (I + dt*H[rho]) u^{n+1} = u^n,   H = -D2(hbar^2/2m) + diag(vloc),
       ! solved by the LAPACK tridiagonal routine dgtsv. The local potential
       ! vloc[rho] is updated by a few inner fixed-point sweeps (re-evaluated at
       ! the provisional new density) so the frozen-operator solve stays
       ! consistent at large dt. Unconditionally stable in the kinetic part, so
       ! dt can be orders of magnitude larger than the explicit-potential
       ! ceiling, reducing the iteration count from thousands to tens.
       call backward_euler_step_l0(nr, r, dr, dt, vloc, phi_in, phi_out)
       block
         integer :: be_it
         real(dp), allocatable :: vloc_new(:), rho_new(:)
         allocate(vloc_new(nr), rho_new(nr))
         do be_it = 1, 6
            call phi_to_density(nr, phi_out, rho_new)
            call compute_ite_vloc_constkf(nr, r, dr, w, rho_new, nmu, nr_tab, rtab, ftab, vloc_new)
            vloc_new(:) = 0.5_dp*(vloc_new(:) + vloc(:))
            call backward_euler_step_l0(nr, r, dr, dt, vloc_new, phi_in, phi_out)
         end do
         deallocate(vloc_new, rho_new)
       end block

    case ('ITE_2AF')
       call potential_step(nr, dt, 0.5_dp, vloc, phi_in, tmp)
       call kinetic_step_l0(nr, r, dt, 1.0_dp, tmp, tmp2)
       call potential_fpi(nr, r, dr, w, tmp2, phi_in, npart_species, dt, 0.5_dp, &
            fpi_max_iter, fpi_tol, nmu, nr_tab, rtab, ftab, phi_out)

    case ('ITE_4A00')
       ! Chin-Krotscheck-like fourth-order step, radial prototype.
       ! The gradient correction follows the uploaded split_time_module structure.
       call potential_step(nr, dt, 1.0_dp/6.0_dp, vloc, phi_in, tmp)
       call kinetic_step_l0(nr, r, dt, 0.5_dp, tmp, tmp2)
       call first_derivative_uniform(nr, dr, vloc, vgrad)
       ! [V,[T,V]] = (hbar^2/m) |grad V|^2 ; the hbar^2/m factor was missing in
       ! the unit-free (hbar=m=1) prototype and is essential in nuclear units.
       vtilde(:) = vloc(:) + dt*dt * (2.0_dp*hbar2_over_2m) * vgrad(:)*vgrad(:) / 48.0_dp
       call potential_step(nr, dt, 2.0_dp/3.0_dp, vtilde, tmp2, tmp)
       call kinetic_step_l0(nr, r, dt, 0.5_dp, tmp, tmp2)
       call potential_step(nr, dt, 1.0_dp/6.0_dp, vloc, tmp2, phi_out)

    case ('ITE_TI')
       ! Takahashi-Imada: Strang splitting with the trace-level fourth-order
       ! corrected potential  Vt = V + (dt^2/24) [V,[T,V]]
       !                        = V + (dt^2/24)(hbar^2/m) |grad V|^2.
       ! All coefficients positive (admissible in imaginary time); one
       ! convolution, two kinetic half-steps and one potential step per
       ! iteration -- same cost as ITE_2A0.
       call first_derivative_uniform(nr, dr, vloc, vgrad)
       vtilde(:) = vloc(:) + dt*dt * (2.0_dp*hbar2_over_2m) * vgrad(:)*vgrad(:) / 24.0_dp
       call kinetic_step_l0(nr, r, dt, 0.5_dp, phi_in, tmp)
       call potential_step(nr, dt, 1.0_dp, vtilde, tmp, tmp2)
       call kinetic_step_l0(nr, r, dt, 0.5_dp, tmp2, phi_out)

    case default
       write(*,*) 'WARNING: unknown ITE method; using ITE_2A0: ', trim(method)
       call potential_step(nr, dt, 0.5_dp, vloc, phi_in, tmp)
       call kinetic_step_l0(nr, r, dt, 1.0_dp, tmp, tmp2)
       call potential_step(nr, dt, 0.5_dp, vloc, tmp2, phi_out)
    end select

    call normalize_phi(nr, r, w, npart_species, phi_out)
    deallocate(rho, vloc, tmp, tmp2, vgrad, vtilde)
  end subroutine ite_step_constkf


  subroutine compute_ite_vloc_constkf(nr, r, dr, w, rho_q, nmu, nr_tab, rtab, ftab, v_eff_ite)
    integer, intent(in) :: nr, nmu, nr_tab
    real(dp), intent(in) :: r(nr), dr, w(nr), rho_q(nr)
    real(dp), intent(in) :: rtab(nr_tab), ftab(nr_tab)
    real(dp), intent(out) :: v_eff_ite(nr)

    real(dp), allocatable :: conv(:), dconv(:), rho_tot(:), v_skp(:), phi(:), lap_phi(:)
    real(dp), allocatable :: v_vw(:), v_kernel(:), v_corr(:), v_sl(:), dTdrho(:), dEdrho(:)
    real(dp), allocatable :: hphi(:), res_phi(:)
    real(dp) :: mu_est, residual_norm

    allocate(conv(nr), dconv(nr), rho_tot(nr), v_skp(nr), phi(nr), lap_phi(nr))
    allocate(v_vw(nr), v_kernel(nr), v_corr(nr), v_sl(nr), dTdrho(nr), dEdrho(nr))
    allocate(hphi(nr), res_phi(nr))

    call compute_kernel_conv(nr, r, w, rho_q, nmu, nr_tab, rtab, ftab, conv, dconv)
    rho_tot(:) = 2.0_dp*rho_q(:)
    call compute_skp_potential(nr, r, dr, rho_tot, v_skp)
    call compute_constkf_derivative_arrays(nr, r, dr, w, rho_q, conv, dconv, v_skp, &
         phi, lap_phi, v_vw, v_kernel, v_corr, v_sl, dTdrho, dEdrho, v_eff_ite, hphi, &
         mu_est, res_phi, residual_norm)

    deallocate(conv, rho_tot, v_skp, phi, lap_phi)
    deallocate(v_vw, v_kernel, v_corr, v_sl, dTdrho, dEdrho, hphi, res_phi)
  end subroutine compute_ite_vloc_constkf


  subroutine evaluate_ite_energies(nr, r, dr, w, rho, nmu, nr_tab, rtab, ftab, &
       T_vw, T_ker, T_corr, T_sl, E_skp, E_sym, mu_est, residual_norm)
    ! Total-energy and Euler-residual monitor for the ITE:
    !   E_sym = 2 (T_vW + T_kernel + T_corr + T_SL) + E_SkP[2 rho].
    integer, intent(in) :: nr, nmu, nr_tab
    real(dp), intent(in) :: r(nr), dr, w(nr), rho(nr)
    real(dp), intent(in) :: rtab(nr_tab), ftab(nr_tab)
    real(dp), intent(out) :: T_vw, T_ker, T_corr, T_sl, E_skp, E_sym, mu_est, residual_norm

    real(dp), allocatable :: conv(:), dconv(:), tau_tf(:), tau_vw(:), edens(:), rho_tot(:), v_skp(:)
    real(dp), allocatable :: phi(:), lap_phi(:), v_vw(:), v_kernel(:), v_corr(:), v_sl(:)
    real(dp), allocatable :: dTdrho(:), dEdrho(:), v_eff(:), hphi(:), res_phi(:)
    integer :: i

    allocate(conv(nr), dconv(nr), tau_tf(nr), tau_vw(nr), edens(nr), rho_tot(nr), v_skp(nr))
    allocate(phi(nr), lap_phi(nr), v_vw(nr), v_kernel(nr), v_corr(nr), v_sl(nr))
    allocate(dTdrho(nr), dEdrho(nr), v_eff(nr), hphi(nr), res_phi(nr))

    call compute_kernel_conv(nr, r, w, rho, nmu, nr_tab, rtab, ftab, conv, dconv)
    call compute_tf_vw_tau(nr, r, dr, rho, tau_tf, tau_vw)

    T_vw = hbar2_over_2m * integrate_r3(nr, r, w, tau_vw)

    do i = 1, nr
       edens(i) = c_NL * max(rho(i),0.0_dp)**(5.0_dp/6.0_dp) * conv(i)
    end do
    T_ker = hbar2_over_2m * integrate_r3(nr, r, w, edens)

    do i = 1, nr
       edens(i) = -c_CORR * max(rho(i),0.0_dp)**(5.0_dp/3.0_dp)
    end do
    T_corr = hbar2_over_2m * integrate_r3(nr, r, w, edens)

    ! Spatial-gradient regularization energy, consistent with the potential
    ! used by the SCF/ITE/diagonalization: E_grad = s2 (hbar^2/2m) int rho^(1/3)|grad rho|^2.
    if (abs(s2_grad) > 0.0_dp) then
       block
         real(dp), allocatable :: rp(:)
         allocate(rp(nr))
         call first_derivative_uniform(nr, dr, rho, rp)
         do i = 1, nr
            edens(i) = s2_grad * max(rho(i),rho_floor)**(1.0_dp/3.0_dp) * rp(i)**2
         end do
         T_corr = T_corr + hbar2_over_2m * integrate_r3(nr, r, w, edens)
         deallocate(rp)
       end block
    end if

    T_sl = sl_energy(nr, r, dr, w, rho)

    rho_tot(:) = 2.0_dp*rho(:)
    E_skp = compute_skp_energy(nr, r, dr, w, rho_tot)
    E_sym = 2.0_dp*(T_vw + T_ker + T_corr + T_sl) + E_skp

    call compute_skp_potential(nr, r, dr, rho_tot, v_skp)
    call compute_constkf_derivative_arrays(nr, r, dr, w, rho, conv, dconv, v_skp, &
         phi, lap_phi, v_vw, v_kernel, v_corr, v_sl, dTdrho, dEdrho, v_eff, hphi, &
         mu_est, res_phi, residual_norm)

    deallocate(conv, tau_tf, tau_vw, edens, rho_tot, v_skp)
    deallocate(phi, lap_phi, v_vw, v_kernel, v_corr, v_sl)
    deallocate(dTdrho, dEdrho, v_eff, hphi, res_phi)
  end subroutine evaluate_ite_energies


  subroutine init_radial_kinetic_propagator(nr, dr)
    integer, intent(in) :: nr
    real(dp), intent(in) :: dr
    real(dp), allocatable :: diag(:), offdiag(:)
    integer :: i

    if (kinetic_ready) then
       if (kinetic_nr_saved == nr .and. abs(kinetic_dr_saved-dr) < 1.0e-14_dp) return
    end if

    if (allocated(kinetic_eval)) deallocate(kinetic_eval)
    if (allocated(kinetic_evec)) deallocate(kinetic_evec)

    allocate(diag(nr), offdiag(nr-1), kinetic_eval(nr), kinetic_evec(nr,nr))

    do i = 1, nr
       diag(i) = 2.0_dp*hbar2_over_2m/(dr*dr)
    end do
    do i = 1, nr-1
       offdiag(i) = -hbar2_over_2m/(dr*dr)
    end do

    call diagonalize_tridiagonal(nr, diag, offdiag, kinetic_eval, kinetic_evec)

    kinetic_ready = .true.
    kinetic_nr_saved = nr
    kinetic_dr_saved = dr
    deallocate(diag, offdiag)
  end subroutine init_radial_kinetic_propagator


  subroutine kinetic_step_l0(nr, r, dt, coeff, phi_in, phi_out)
    integer, intent(in) :: nr
    real(dp), intent(in) :: r(nr), dt, coeff, phi_in(nr)
    real(dp), intent(out) :: phi_out(nr)
    real(dp), allocatable :: u(:), c(:), uout(:)
    integer :: i, j
    real(dp) :: factor

    if (.not. kinetic_ready) stop 'ERROR: kinetic propagator is not initialized.'

    allocate(u(nr), c(nr), uout(nr))
    do i = 1, nr
       u(i) = r(i)*phi_in(i)
    end do

    c(:) = 0.0_dp
    do j = 1, nr
       do i = 1, nr
          c(j) = c(j) + kinetic_evec(i,j)*u(i)
       end do
       factor = exp(max(-80.0_dp, min(80.0_dp, -coeff*dt*kinetic_eval(j))))
       c(j) = c(j)*factor
    end do

    uout(:) = 0.0_dp
    do j = 1, nr
       do i = 1, nr
          uout(i) = uout(i) + kinetic_evec(i,j)*c(j)
       end do
    end do

    do i = 1, nr
       phi_out(i) = uout(i)/max(r(i), 1.0e-12_dp)
    end do
    deallocate(u, c, uout)
  end subroutine kinetic_step_l0


  subroutine potential_implicit_step(nr, r, dr, w, dt, phi_in, nmu, nr_tab, rtab, ftab, phi_out)
    ! Implicit (backward-Euler) potential-only step for the diagonal local
    ! potential:  phi_out = phi_in / (1 + dt*vloc[rho_out]),  with rho_out
    ! = phi_out^2 resolved by a few fixed-point sweeps. The kinetic operator
    ! is NOT touched here (handled by the surrounding kinetic half-steps in
    ! the Strang composition). Unconditionally stable for vloc >= 0 and, with
    ! the clamp below, robust for the mildly negative local potential too.
    integer, intent(in) :: nr, nmu, nr_tab
    real(dp), intent(in) :: r(nr), dr, w(nr), dt, phi_in(nr)
    real(dp), intent(in) :: rtab(nr_tab), ftab(nr_tab)
    real(dp), intent(out) :: phi_out(nr)
    real(dp), allocatable :: rho(:), vloc(:), denom(:)
    integer :: it, i

    allocate(rho(nr), vloc(nr), denom(nr))
    phi_out(:) = phi_in(:)
    do it = 1, 6
       call phi_to_density(nr, phi_out, rho)
       call compute_ite_vloc_constkf(nr, r, dr, w, rho, nmu, nr_tab, rtab, ftab, vloc)
       do i = 1, nr
          ! backward-Euler denominator; guard against 1+dt*v <= 0 for strongly
          ! attractive v by flooring at a small positive value.
          denom(i) = max(1.0_dp + dt*vloc(i), 1.0e-3_dp)
          phi_out(i) = phi_in(i)/denom(i)
       end do
    end do
    deallocate(rho, vloc, denom)
  end subroutine potential_implicit_step


  subroutine backward_euler_step_l0(nr, r, dr, dt, vloc, phi_in, phi_out)
    ! One fully-implicit (backward-Euler) imaginary-time step on u = r*phi:
    !   (I + dt*H) u_out = u_in,  H = -(hbar^2/2m) d^2/dr^2 + diag(vloc),
    ! with Dirichlet ends, solved by LAPACK dgtsv. H is frozen at the incoming
    ! density (vloc supplied by the caller); the step is unconditionally stable.
    integer, intent(in) :: nr
    real(dp), intent(in) :: r(nr), dr, dt, vloc(nr), phi_in(nr)
    real(dp), intent(out) :: phi_out(nr)
    real(dp), allocatable :: dl(:), dd(:), du(:), bb(:)
    real(dp) :: k
    integer :: i, info

    allocate(dl(nr-1), dd(nr), du(nr-1), bb(nr))
    k = hbar2_over_2m/(dr*dr)
    ! (I + dt*H): diagonal 1 + dt*(2k + vloc); off-diagonals -dt*k
    do i = 1, nr
       dd(i) = 1.0_dp + dt*(2.0_dp*k + vloc(i))
       bb(i) = r(i)*phi_in(i)
    end do
    do i = 1, nr-1
       dl(i) = -dt*k
       du(i) = -dt*k
    end do
    call dgtsv(nr, 1, dl, dd, du, bb, nr, info)
    if (info /= 0) then
       write(*,*) 'ERROR: backward_euler_step_l0 dgtsv info=', info
       stop
    end if
    do i = 1, nr
       phi_out(i) = bb(i)/max(r(i), 1.0e-12_dp)
    end do
    deallocate(dl, dd, du, bb)
  end subroutine backward_euler_step_l0


  subroutine potential_step(nr, dt, coeff, vloc, phi_in, phi_out)
    integer, intent(in) :: nr
    real(dp), intent(in) :: dt, coeff, vloc(nr), phi_in(nr)
    real(dp), intent(out) :: phi_out(nr)
    integer :: i
    real(dp) :: arg

    do i = 1, nr
       arg = -coeff*dt*vloc(i)
       arg = max(-80.0_dp, min(80.0_dp, arg))
       phi_out(i) = exp(arg)*phi_in(i)
    end do
  end subroutine potential_step


  subroutine potential_fpi(nr, r, dr, w, phi_base, phi_init, npart_species, dt, coeff, &
       max_fpi, tol_fpi, nmu, nr_tab, rtab, ftab, phi_out)
    integer, intent(in) :: nr, max_fpi, nmu, nr_tab
    real(dp), intent(in) :: r(nr), dr, w(nr), phi_base(nr), phi_init(nr), npart_species, dt, coeff, tol_fpi
    real(dp), intent(in) :: rtab(nr_tab), ftab(nr_tab)
    real(dp), intent(out) :: phi_out(nr)

    real(dp), allocatable :: phi0(:), phi1(:), rho(:), vloc(:)
    integer :: iter
    real(dp) :: err

    allocate(phi0(nr), phi1(nr), rho(nr), vloc(nr))
    phi0(:) = phi_init(:)
    call normalize_phi(nr, r, w, npart_species, phi0)

    do iter = 1, max_fpi
       call phi_to_density(nr, phi0, rho)
       call compute_ite_vloc_constkf(nr, r, dr, w, rho, nmu, nr_tab, rtab, ftab, vloc)
       call potential_step(nr, dt, coeff, vloc, phi_base, phi1)
       call normalize_phi(nr, r, w, npart_species, phi1)
       err = phi_l2_difference(nr, r, w, phi0, phi1)
       phi0(:) = phi1(:)
       if (err < tol_fpi) exit
    end do
    phi_out(:) = phi0(:)
    deallocate(phi0, phi1, rho, vloc)
  end subroutine potential_fpi


  subroutine normalize_phi(nr, r, w, npart_species, phi)
    integer, intent(in) :: nr
    real(dp), intent(in) :: r(nr), w(nr), npart_species
    real(dp), intent(inout) :: phi(nr)
    real(dp) :: rho_tmp(nr), norm
    integer :: i

    do i = 1, nr
       rho_tmp(i) = max(phi(i)*phi(i), rho_floor)
    end do
    norm = integrate_r3(nr, r, w, rho_tmp)
    if (norm <= 0.0_dp) stop 'ERROR: phi normalization failed.'
    phi(:) = phi(:)*sqrt(npart_species/norm)
  end subroutine normalize_phi


  real(dp) function phi_l2_difference(nr, r, w, phi_a, phi_b)
    integer, intent(in) :: nr
    real(dp), intent(in) :: r(nr), w(nr), phi_a(nr), phi_b(nr)
    integer :: i
    real(dp) :: tmp

    phi_l2_difference = 0.0_dp
    do i = 1, nr
       tmp = phi_b(i) - phi_a(i)
       phi_l2_difference = phi_l2_difference + 4.0_dp*pi*r(i)**2 * tmp*tmp * w(i)
    end do
    phi_l2_difference = sqrt(max(phi_l2_difference, 0.0_dp))
  end function phi_l2_difference


  subroutine write_density_output(filename, nr, r, rho)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nr
    real(dp), intent(in) :: r(nr), rho(nr)
    integer :: i, unit

    unit = 82
    open(unit=unit, file=trim(filename), status='replace', action='write')
    write(unit,'(A)') '# r[fm] rho_q_ITE[fm^-3]'
    do i = 1, nr
       write(unit,'(2ES24.15)') r(i), rho(i)
    end do
    close(unit)
  end subroutine write_density_output

end module radial_ite_split_module

module diag_scf_module
  ! Rearranged-diagonalization SCF for the OF Euler equation, following
  ! Imoto, Imada & Oshiyama, Phys. Rev. Research 3, 033198 (2021), Eq. (16):
  !   T = lambda T_vW + T_r,
  !   [ -C d^2/dr^2 + (1/lambda)( (1-lambda) v_vW[rho_n] + v_rest[rho_n] ) ] u
  !        = (mu/lambda) u,   u = r sqrt(rho), l=0,
  ! solved by tridiagonal diagonalization (LAPACK dstev), lowest eigenstate,
  ! then linear density mixing rho <- (1-alpha) rho + alpha rho_new.
  ! lambda = 1 is the naive vW-partitioned scheme (known to converge poorly);
  ! lambda >> 1 damps the potential feedback.  The fixed point is
  ! lambda-independent.
  use precision_module
  use constants_module
  use lindhard_kernel_module, only: build_f1_table, check_f1_norm
  use lapack_tridiag_module, only: dstev
  use radial_ite_split_module, only: compute_ite_vloc_constkf, evaluate_ite_energies
  implicit none
contains

  subroutine run_diag_scf(nr, r, dr, w, rho_init, xNpart, kf_const, qmax, nq, nmu, nr_tab, &
       lam, alpha_mix, max_iter, tol, history_file, density_file, summary_file, &
       scf_mode, inv_shift)
    ! Rearranged-diagonalization SCF (Imoto-Imada-Oshiyama). The frozen local
    ! potential is built by compute_ite_vloc_constkf, which internally calls the
    ! kernel dispatcher compute_kernel_conv; that dispatcher is dd-aware, so this
    ! solver handles BOTH the constant-kF benchmark AND the full
    ! density-dependent kernel (with the rearrangement term) whenever
    ! kernel_mode='dd' has been set. The dd kernel typically needs stronger
    ! mixing (alpha ~ 0.1-0.5) than constant-kF (alpha = 1) at equal gap.
    ! scf_mode = 'eig'   : exact lowest eigenstate of the frozen H (dstev), i.e.
    !                      fully converged inner loop per outer step.
    ! scf_mode = 'invit' : Henning-type BASIC INVERSE ITERATION -- one solve of
    !                      (H_lambda[rho_n] + s) w = u_n per outer step
    !                      (positivity shift s = inv_shift, required in the
    !                      self-bound nuclear case where lambda_1 = mu/lambda < 0
    !                      violates Henning's coercivity assumptions).
    integer, intent(in) :: nr, nq, nmu, nr_tab, max_iter
    real(dp), intent(in) :: r(nr), dr, w(nr), rho_init(nr), xNpart, kf_const, qmax
    real(dp), intent(in) :: lam, alpha_mix, tol
    character(len=*), intent(in) :: history_file, density_file, summary_file
    character(len=*), intent(in) :: scf_mode
    real(dp), intent(in) :: inv_shift

    real(dp), allocatable :: rtab(:), ftab(:), vol(:)
    real(dp), allocatable :: rho(:), rho_new(:), u(:), vrest(:), vvw(:), veff(:)
    real(dp), allocatable :: dmat(:), emat(:), zmat(:,:), work(:)
    real(dp) :: rtab_max, cC, s2, resid, eps0, eps1, mu_est, Ecur, dE
    integer :: i, it, info, unit
    logical :: converged

    allocate(rtab(nr_tab), ftab(nr_tab), vol(nr))
    allocate(rho(nr), rho_new(nr), u(nr), vrest(nr), vvw(nr), veff(nr))
    allocate(dmat(nr), emat(nr-1), zmat(nr,nr), work(max(1,2*nr-2)))

    cC = hbar2_over_2m
    rtab_max = 2.0_dp * r(nr)
    call build_f1_table(kf_const, qmax, nq, nr_tab, rtab_max, rtab, ftab)
    call check_f1_norm(nr_tab, rtab, ftab)

    do i = 1, nr
       vol(i) = 4.0_dp*pi*r(i)*r(i)*w(i)
    end do
    rho(:) = max(rho_init(:), 0.0_dp)
    u(:) = r(:)*sqrt(rho(:))

    open(newunit=unit, file=history_file, status='replace', action='write')
    write(unit,'(A)') '# iter  E_sym  mu_est  resid_rho  eps0  gap(eps1-eps0)'

    converged = .false.
    Ecur = 0.0_dp
    do it = 1, max_iter
       ! frozen potentials at the current density
       call compute_ite_vloc_constkf(nr, r, dr, w, rho, nmu, nr_tab, rtab, ftab, vrest)
       do i = 1, nr
          if (rho(i) > 1.0e-12_dp .and. abs(u(i)) > 1.0e-14_dp) then
             if (i == 1) then
                vvw(i) = -cC*(u(2) - 2.0_dp*u(1))/(dr*dr)/u(1)
             else if (i == nr) then
                vvw(i) = -cC*(-2.0_dp*u(nr) + u(nr-1))/(dr*dr)/u(nr)
             else
                vvw(i) = -cC*(u(i+1) - 2.0_dp*u(i) + u(i-1))/(dr*dr)/u(i)
             end if
             vvw(i) = max(min(vvw(i), 500.0_dp), -500.0_dp)
          else
             vvw(i) = 0.0_dp
          end if
          veff(i) = ((1.0_dp - lam)*vvw(i) + vrest(i))/lam
       end do

       ! tridiagonal H = -C D2 + diag(veff) on u, Dirichlet BCs
       ! (eigenvalues monitored in both modes)
       do i = 1, nr
          dmat(i) = 2.0_dp*cC/(dr*dr) + veff(i)
       end do
       emat(:) = -cC/(dr*dr)
       if (trim(scf_mode) == 'invit') then
          block
            real(dp), allocatable :: dl(:), dd(:), du(:), bb(:,:), dsave(:), esave(:)
            integer :: infi
            allocate(dl(nr-1), dd(nr), du(nr-1), bb(nr,1), dsave(nr), esave(nr-1))
            dsave = dmat; esave = emat
            dd(:) = dmat(:) + inv_shift
            dl(:) = emat(:)
            du(:) = emat(:)
            bb(:,1) = u(:)
            call dgtsv(nr, 1, dl, dd, du, bb, nr, infi)
            if (infi /= 0) then
               write(*,*) 'dgtsv failed, info=', infi
               deallocate(dl, dd, du, bb, dsave, esave)
               exit
            end if
            u(:) = bb(:,1)
            ! monitor eigenvalues of the same frozen H
            call dstev('N', nr, dsave, esave, zmat, nr, work, info)
            eps0 = dsave(1)
            eps1 = dsave(2)
            deallocate(dl, dd, du, bb, dsave, esave)
          end block
       else
          call dstev('V', nr, dmat, emat, zmat, nr, work, info)
          if (info /= 0) then
             write(*,*) 'dstev failed, info=', info
             exit
          end if
          eps0 = dmat(1)
          eps1 = dmat(2)
          u(:) = zmat(:,1)
       end if
       if (sum(u) < 0.0_dp) u = -u
       do i = 1, nr
          rho_new(i) = (u(i)/r(i))**2
       end do
       ! normalize to particle number
       s2 = 0.0_dp
       do i = 1, nr
          s2 = s2 + vol(i)*rho_new(i)
       end do
       rho_new(:) = rho_new(:)*xNpart/max(s2,1.0e-300_dp)

       resid = 0.0_dp
       do i = 1, nr
          resid = resid + vol(i)*(rho_new(i)-rho(i))**2
       end do
       resid = sqrt(resid)

       rho(:) = (1.0_dp - alpha_mix)*rho(:) + alpha_mix*rho_new(:)
       u(:) = r(:)*sqrt(max(rho(:),0.0_dp))

       block
         real(dp) :: Tvw_m, Tker_m, Tcorr_m, Tsl_m, Eskp_m, Esym_m, mu_m, res_m
         call evaluate_ite_energies(nr, r, dr, w, rho, nmu, nr_tab, rtab, ftab, &
              Tvw_m, Tker_m, Tcorr_m, Tsl_m, Eskp_m, Esym_m, mu_m, res_m)
         dE = Esym_m - Ecur
         Ecur = Esym_m
       end block
       mu_est = lam*eps0

       write(unit,'(I8,5ES20.10)') it, Ecur, mu_est, resid, eps0, eps1-eps0
       if (resid < tol .and. it > 3) then
          converged = .true.
          exit
       end if
    end do
    close(unit)

    open(newunit=unit, file=density_file, status='replace', action='write')
    write(unit,'(A)') '# r[fm] rho_q_DIAG[fm^-3]'
    do i = 1, nr
       write(unit,'(2ES24.15)') r(i), rho(i)
    end do
    close(unit)

    open(newunit=unit, file=summary_file, status='replace', action='write')
    write(unit,'(A,ES24.15)') 'E_total_NeqZ_MeV        ', Ecur
    write(unit,'(A,ES24.15)') 'mu_est_MeV              ', mu_est
    write(unit,'(A,ES24.15)') 'resid_rho               ', resid
    write(unit,'(A,F10.4)')   'lambda                  ', lam
    write(unit,'(A,F10.4)')   'alpha_mix               ', alpha_mix
    write(unit,'(A,L2)')      'converged               ', converged
    close(unit)
    write(*,'(A,ES16.8,A,ES11.3,A,L2)') ' Diag-SCF done: E = ', Ecur, '  resid = ', resid, '  conv = ', converged

    deallocate(rtab, ftab, vol, rho, rho_new, u, vrest, vvw, veff, dmat, emat, zmat, work)
  end subroutine run_diag_scf

end module diag_scf_module


program ks_skp_plus_nlf_constkf_v8
  use precision_module
  use constants_module
  use radial_grid_module
  use ks_scf_module
  use density_file_module
  use nlf_postprocess_module
  use nl_ofdft_derivative_module
  use radial_initial_density_module
  use radial_ite_split_module
  use semilocal_correction_module
  use dd_teacher_module
  use derivative_checker_module
  use direct_minimization_module
  use diag_scf_module
  implicit none

  ! Main switches.
  logical :: run_ks, run_nlf, run_derivative, run_ite, nlf_use_ks_density

  ! KS parameters.
  integer :: A, npart_q, nr, lmax, nstates_per_l, max_scf
  real(dp) :: rbox, scf_tol, mix_alpha, init_radius, init_diffuseness

  ! NLF parameters.
  real(dp) :: kf_const, qmax, density_scale
  integer :: nq, nmu, nr_tab

  ! ITE parameters.
  integer :: ite_max_iter, ite_fpi_max_iter, ite_monitor_every
  real(dp) :: ite_dt, ite_tol, ite_fpi_tol, ite_init_r0, ite_init_a
  logical :: ite_dt_adapt
  real(dp) :: ite_dt_max, ite_dt_grow, ite_dt_shrink
  character(len=32) :: ite_method, ite_init_mode

  ! Reference-data generator, semilocal-fit, and derivative-checker parameters.
  ! The public namelist keeps the legacy run_teacher/teacher_* names for
  ! compatibility with existing input files.
  logical :: run_teacher, run_deriv_check, use_semilocal, sl_ok
  real(dp) :: teacher_weight_alpha, teacher_xqmax
  integer :: teacher_ny_gtab
  real(dp) :: check_eta_r0, check_eta_sigma
  character(len=256) :: sl_coeff_file, teacher_file, teacher_dataset_file, deriv_check_file
  character(len=8) :: kernel_mode
  real(dp) :: dd_kf_scale_in, dd_xqmax_in
  integer :: dd_ny_base_in
  real(dp) :: kernel_lambda
  logical :: run_minimize
  integer :: min_max_iter
  real(dp) :: min_gtol
  character(len=256) :: min_history_file, min_density_file, min_summary_file
  character(len=256) :: extra_kernel_file
  character(len=256) :: extra_g_file
  real(dp) :: c_corr_shift
  real(dp) :: s2_grad_in
  logical :: run_diagscf
  real(dp) :: diag_lambda, diag_alpha, diag_tol
  integer :: diag_max_iter
  character(len=256) :: diag_history_file, diag_density_file, diag_summary_file
  character(len=8) :: diag_mode
  real(dp) :: diag_shift

  ! File names.
  character(len=256) :: density_species_file, density_total_file
  character(len=256) :: levels_file, history_file, ks_nlf_file
  character(len=256) :: external_density_file
  character(len=256) :: nlf_file, kernel_file
  character(len=256) :: derivative_file, derivative_energy_file, derivative_kernel_file
  character(len=256) :: ite_initial_density_file, ite_density_file, ite_history_file
  character(len=256) :: ite_derivative_file, ite_summary_file, ite_kernel_file

  ! Arrays.
  integer :: ios, nr_ext
  real(dp) :: dr, dr_ext
  real(dp), allocatable :: r(:), w(:), rho_q(:), rho_tot(:)
  real(dp), allocatable :: r_ext(:), rho_ext(:), w_ext(:)
  real(dp), allocatable :: rho_ite_init(:)

  namelist /input/ run_ks, run_nlf, nlf_use_ks_density, &
       A, npart_q, nr, rbox, lmax, nstates_per_l, max_scf, scf_tol, mix_alpha, &
       init_radius, init_diffuseness, &
       density_species_file, density_total_file, levels_file, history_file, ks_nlf_file, &
       external_density_file, density_scale, &
       kf_const, qmax, nq, nmu, nr_tab, nlf_file, kernel_file, &
       run_derivative, derivative_file, derivative_energy_file, derivative_kernel_file, &
       run_ite, ite_method, ite_init_mode, ite_max_iter, ite_dt, ite_tol, &
       ite_fpi_max_iter, ite_fpi_tol, ite_init_r0, ite_init_a, ite_initial_density_file, &
       ite_density_file, ite_history_file, ite_derivative_file, ite_summary_file, ite_kernel_file, &
       ite_monitor_every, &
       ite_dt_adapt, ite_dt_max, ite_dt_grow, ite_dt_shrink, &
       run_teacher, teacher_weight_alpha, teacher_xqmax, teacher_ny_gtab, &
       teacher_file, teacher_dataset_file, &
       use_semilocal, sl_coeff_file, &
       run_deriv_check, check_eta_r0, check_eta_sigma, deriv_check_file, &
       kernel_mode, dd_kf_scale_in, dd_xqmax_in, dd_ny_base_in, kernel_lambda, &
       run_minimize, min_max_iter, min_gtol, min_history_file, min_density_file, min_summary_file, &
       extra_kernel_file, c_corr_shift, extra_g_file, s2_grad_in, &
       run_diagscf, diag_lambda, diag_alpha, diag_max_iter, diag_tol, &
       diag_history_file, diag_density_file, diag_summary_file, &
       diag_mode, diag_shift

  ! Defaults: 16O test.
  run_ks = .true.
  run_nlf = .true.
  run_derivative = .true.
  run_ite = .false.
  nlf_use_ks_density = .true.

  A  = 16
  npart_q = 8
  nr = 500
  rbox = 20.0_dp
  lmax = 5
  nstates_per_l = 4
  max_scf = 200
  scf_tol = 1.0e-7_dp
  mix_alpha = 0.20_dp
  init_radius = -1.0_dp
  init_diffuseness = 0.55_dp

  density_species_file = 'ks_density_species.dat'
  density_total_file   = 'ks_density_total.dat'
  levels_file          = 'ks_levels.dat'
  history_file         = 'ks_scf_history.dat'
  ks_nlf_file          = 'ks_nlf.dat'

  external_density_file = 'ks_density_species.dat'
  density_scale = 1.0_dp

  kf_const = 0.80_dp
  qmax = 60.0_dp
  nq = 6000
  nmu = 32
  nr_tab = 4001
  nlf_file = 'nlf_constkf.dat'
  kernel_file = 'kernel_f1_constkf.dat'
  derivative_file = 'ofdft_derivative_constkf.dat'
  derivative_energy_file = 'ofdft_derivative_summary.dat'
  derivative_kernel_file = 'kernel_f1_constkf_derivative.dat'

  ite_method = 'ITE_2A0'
  ite_init_mode = 'ks'
  ite_max_iter = 20
  ite_dt = 1.0e-4_dp
  ite_tol = 1.0e-8_dp
  ite_fpi_max_iter = 50
  ite_fpi_tol = 1.0e-12_dp
  ite_init_r0 = 1.20_dp
  ite_init_a = 0.55_dp
  ite_monitor_every = 1
  ite_dt_adapt = .false.
  ite_dt_max = 5.0e-5_dp
  ite_dt_grow = 1.1_dp
  ite_dt_shrink = 0.5_dp
  run_teacher = .false.
  teacher_weight_alpha = 1.5_dp
  teacher_xqmax = 80.0_dp
  teacher_ny_gtab = 8001
  teacher_file = 'teacher_ddkf.dat'
  teacher_dataset_file = 'fit_dataset.dat'
  use_semilocal = .false.
  sl_coeff_file = 'fit_coefficients.dat'
  run_deriv_check = .false.
  check_eta_r0 = 2.5_dp
  check_eta_sigma = 1.0_dp
  deriv_check_file = 'derivative_check.dat'
  kernel_mode = 'const'
  dd_kf_scale_in = 1.0_dp
  dd_xqmax_in = 80.0_dp
  dd_ny_base_in = 8001
  kernel_lambda = 1.0_dp
  run_minimize = .false.
  min_max_iter = 5000
  min_gtol = 1.0d-8
  min_history_file = 'min_history.dat'
  min_density_file = 'min_density.dat'
  min_summary_file = 'min_summary.dat'
  extra_kernel_file = ''
  extra_g_file = ''
  c_corr_shift = 0.0d0
  s2_grad_in = 0.0d0
  run_diagscf = .false.
  diag_lambda = 10.0d0
  diag_alpha = 0.3d0
  diag_max_iter = 500
  diag_tol = 1.0d-10
  diag_history_file = 'diag_history.dat'
  diag_density_file = 'diag_density.dat'
  diag_summary_file = 'diag_summary.dat'
  diag_mode = 'eig'
  diag_shift = 30.0d0

  ite_initial_density_file = 'ks_density_species.dat'
  ite_density_file = 'ite_density_constkf.dat'
  ite_history_file = 'ite_history_constkf.dat'
  ite_derivative_file = 'ite_derivative_constkf.dat'
  ite_summary_file = 'ite_derivative_summary.dat'
  ite_kernel_file = 'ite_kernel_f1_constkf.dat'

  open(unit=99, file='input_ks_nlf.dat', status='old', action='read', iostat=ios)
  if (ios == 0) then
     read(99, nml=input)
  ! --- kernel mode ---
  if (trim(adjustl(kernel_mode)) == 'dd') then
     dd_kernel_mode = .true.
     dd_xqmax = dd_xqmax_in
     dd_nq = nq
     dd_ny_base = dd_ny_base_in
     dd_kf_scale = dd_kf_scale_in
     write(*,'(A)') 'Kernel mode: density-dependent kF (variational, exact derivative).'
     write(*,'(A,F10.5)') 'dd kF-scale lambda = ', dd_kf_scale
     write(*,'(A,F8.2,A,I0,A,I0)') 'dd table: xqmax=', dd_xqmax, &
          ' nq=', dd_nq, ' ny_base=', dd_ny_base
  else
     dd_kernel_mode = .false.
     write(*,'(A)') 'Kernel mode: constant kF.'
  end if
  call set_kernel_lambda(kernel_lambda)
  if (len_trim(extra_kernel_file) > 0) then
     call load_extra_kernel(trim(extra_kernel_file))
  end if
  if (len_trim(extra_g_file) > 0) then
     call load_extra_g(trim(extra_g_file))
  end if
  if (abs(c_corr_shift) > 0.0d0) then
     c_CORR = c_CORR + c_corr_shift
     write(*,'(A,ES14.6,A,ES14.6)') ' c_CORR shifted by ', c_corr_shift, ' -> ', c_CORR
  end if
  s2_grad = s2_grad_in
  if (abs(s2_grad) > 0.0d0) then
     write(*,'(A,ES14.6)') ' spatial-gradient regularization s2_grad = ', s2_grad
  end if

     close(99)
  else
     write(*,*) 'input_ks_nlf.dat not found. Using default 16O settings.'
  end if

  if (use_semilocal) then
     call sl_load_coefficients(sl_coeff_file, sl_ok)
     if (.not. sl_ok) then
        write(*,*) 'ERROR: use_semilocal=.true. but the coefficient file could not be loaded.'
        stop
     end if
  end if

  if (run_ks) then
     call run_spherical_ks_skp(A, npart_q, nr, rbox, lmax, nstates_per_l, &
          max_scf, scf_tol, mix_alpha, init_radius, init_diffuseness, &
          density_species_file, density_total_file, levels_file, history_file, ks_nlf_file, &
          r, dr, w, rho_q, rho_tot)
  end if

  if (run_nlf) then
     ! Density source: KS density (if run_ks produced one and nlf_use_ks_density)
     ! or an external density file. Convolution kernel: dd or constant-kF, chosen
     ! by kernel_mode -- INDEPENDENTLY of the density source. This prevents the
     ! dd NLF from silently falling back to the constant-kF postprocessor when the
     ! caller uses the KS density; this is independent of the kernel-mode dispatch.
     if (nlf_use_ks_density) then
        if (.not. allocated(rho_q)) then
           write(*,*) 'ERROR: nlf_use_ks_density=.true. but KS was not run.'
           stop
        end if
        if (trim(adjustl(kernel_mode)) == 'dd') then
           call run_nlf_dd(nr, r, dr, w, rho_q, dd_xqmax_in, nq, dd_ny_base_in, nmu, nlf_file)
        else
           call run_nlf_constkf(nr, r, dr, w, rho_q, kf_const, qmax, nq, nmu, nr_tab, &
                nlf_file, kernel_file)
        end if
     else
        call read_density_2col(external_density_file, density_scale, nr_ext, r_ext, rho_ext)
        allocate(w_ext(nr_ext))
        call make_weights_from_existing_grid(nr_ext, r_ext, w_ext, dr_ext)
        if (trim(adjustl(kernel_mode)) == 'dd') then
           call run_nlf_dd(nr_ext, r_ext, dr_ext, w_ext, rho_ext, dd_xqmax_in, nq, &
                dd_ny_base_in, nmu, nlf_file)
        else
           call run_nlf_constkf(nr_ext, r_ext, dr_ext, w_ext, rho_ext, kf_const, qmax, nq, nmu, nr_tab, &
                nlf_file, kernel_file)
        end if
        deallocate(r_ext, rho_ext, w_ext)
     end if
  end if

  if (run_teacher) then
     if (nlf_use_ks_density) then
        if (.not. allocated(rho_q)) then
           write(*,*) 'ERROR: run_teacher with nlf_use_ks_density=.true. but KS was not run.'
           stop
        end if
        call run_dd_teacher(nr, r, dr, w, rho_q, kf_const, qmax, nq, nmu, nr_tab, &
             teacher_xqmax, teacher_ny_gtab, teacher_weight_alpha, &
             teacher_file, teacher_dataset_file)
     else
        call read_density_2col(external_density_file, density_scale, nr_ext, r_ext, rho_ext)
        allocate(w_ext(nr_ext))
        call make_weights_from_existing_grid(nr_ext, r_ext, w_ext, dr_ext)
        call run_dd_teacher(nr_ext, r_ext, dr_ext, w_ext, rho_ext, kf_const, qmax, nq, nmu, nr_tab, &
             teacher_xqmax, teacher_ny_gtab, teacher_weight_alpha, &
             teacher_file, teacher_dataset_file)
        deallocate(r_ext, rho_ext, w_ext)
     end if
  end if

  if (run_deriv_check) then
     if (nlf_use_ks_density) then
        if (.not. allocated(rho_q)) then
           write(*,*) 'ERROR: run_deriv_check with nlf_use_ks_density=.true. but KS was not run.'
           stop
        end if
        call run_derivative_checker(nr, r, dr, w, rho_q, kf_const, qmax, nq, nmu, nr_tab, &
             check_eta_r0, check_eta_sigma, deriv_check_file)
     else
        call read_density_2col(external_density_file, density_scale, nr_ext, r_ext, rho_ext)
        allocate(w_ext(nr_ext))
        call make_weights_from_existing_grid(nr_ext, r_ext, w_ext, dr_ext)
        call run_derivative_checker(nr_ext, r_ext, dr_ext, w_ext, rho_ext, kf_const, qmax, nq, nmu, nr_tab, &
             check_eta_r0, check_eta_sigma, deriv_check_file)
        deallocate(r_ext, rho_ext, w_ext)
     end if
  end if

  if (run_derivative) then
     if (nlf_use_ks_density) then
        if (.not. allocated(rho_q)) then
           write(*,*) 'ERROR: nlf_use_ks_density=.true. but KS was not run.'
           stop
        end if
        call run_ofdft_derivative_constkf(nr, r, dr, w, rho_q, kf_const, qmax, nq, nmu, nr_tab, &
             derivative_file, derivative_energy_file, derivative_kernel_file)
     else
        call read_density_2col(external_density_file, density_scale, nr_ext, r_ext, rho_ext)
        allocate(w_ext(nr_ext))
        call make_weights_from_existing_grid(nr_ext, r_ext, w_ext, dr_ext)
        call run_ofdft_derivative_constkf(nr_ext, r_ext, dr_ext, w_ext, rho_ext, &
             kf_const, qmax, nq, nmu, nr_tab, derivative_file, derivative_energy_file, &
             derivative_kernel_file)
        deallocate(r_ext, rho_ext, w_ext)
     end if
  end if

  if (run_diagscf) then
     call read_density_2col(external_density_file, density_scale, nr_ext, r_ext, rho_ext)
     allocate(w_ext(nr_ext))
     call make_weights_from_existing_grid(nr_ext, r_ext, w_ext, dr_ext)
     call run_diag_scf(nr_ext, r_ext, dr_ext, w_ext, rho_ext, dble(npart_q), &
          kf_const, qmax, nq, nmu, nr_tab, diag_lambda, diag_alpha, diag_max_iter, diag_tol, &
          diag_history_file, diag_density_file, diag_summary_file, diag_mode, diag_shift)
     deallocate(r_ext, rho_ext, w_ext)
  end if

  if (run_minimize) then
     if (nlf_use_ks_density) then
        if (.not. allocated(rho_q)) then
           write(*,*) 'ERROR: run_minimize with nlf_use_ks_density=.true. but KS was not run.'
           stop
        end if
        call run_direct_minimization(nr, r, dr, w, rho_q, dble(npart_q), kf_const, qmax, nq, nmu, nr_tab, &
             min_max_iter, min_gtol, min_history_file, min_density_file, min_summary_file)
     else
        call read_density_2col(external_density_file, density_scale, nr_ext, r_ext, rho_ext)
        allocate(w_ext(nr_ext))
        call make_weights_from_existing_grid(nr_ext, r_ext, w_ext, dr_ext)
        call run_direct_minimization(nr_ext, r_ext, dr_ext, w_ext, rho_ext, dble(npart_q), &
             kf_const, qmax, nq, nmu, nr_tab, &
             min_max_iter, min_gtol, min_history_file, min_density_file, min_summary_file)
        deallocate(r_ext, rho_ext, w_ext)
     end if
  end if

  if (run_ite) then
     call prepare_and_run_ite()
  end if

  write(*,*) 'All requested tasks finished.'

contains

  subroutine prepare_and_run_ite()
    character(len=32) :: mode

    mode = adjustl(ite_init_mode)

    if (trim(mode) == 'file') then
       call read_density_2col(ite_initial_density_file, density_scale, nr_ext, r_ext, rho_ext)
       allocate(w_ext(nr_ext))
       call make_weights_from_existing_grid(nr_ext, r_ext, w_ext, dr_ext)
       call renormalize_species_density(nr_ext, r_ext, w_ext, real(npart_q,dp), rho_ext)
       call run_ite_constkf_radial(nr_ext, r_ext, dr_ext, w_ext, rho_ext, real(npart_q,dp), &
            ite_method, ite_dt, ite_max_iter, ite_tol, ite_fpi_max_iter, ite_fpi_tol, &
            kf_const, qmax, nq, nmu, nr_tab, ite_monitor_every, ite_density_file, ite_history_file, &
            ite_derivative_file, ite_summary_file, ite_kernel_file, &
            ite_dt_adapt, ite_dt_max, ite_dt_grow, ite_dt_shrink)
       deallocate(r_ext, rho_ext, w_ext)

    else
       if (.not. allocated(r)) then
          allocate(r(nr), w(nr))
          call make_uniform_interior_grid(nr, rbox, r, dr, w)
       end if
       allocate(rho_ite_init(nr))

       if (trim(mode) == 'ks') then
          if (.not. allocated(rho_q)) then
             write(*,*) 'ERROR: ite_init_mode="ks" requires run_ks=.true. in this prototype.'
             stop
          end if
          rho_ite_init(:) = rho_q(:)
       else if (trim(mode) == 'fermi') then
          call init_species_fermi_density(nr, r, w, real(A,dp), real(npart_q,dp), &
               ite_init_r0, ite_init_a, rho_ite_init)
       else
          write(*,*) 'WARNING: unknown ite_init_mode; using Fermi density: ', trim(ite_init_mode)
          call init_species_fermi_density(nr, r, w, real(A,dp), real(npart_q,dp), &
               ite_init_r0, ite_init_a, rho_ite_init)
       end if

       call run_ite_constkf_radial(nr, r, dr, w, rho_ite_init, real(npart_q,dp), &
            ite_method, ite_dt, ite_max_iter, ite_tol, ite_fpi_max_iter, ite_fpi_tol, &
            kf_const, qmax, nq, nmu, nr_tab, ite_monitor_every, ite_density_file, ite_history_file, &
            ite_derivative_file, ite_summary_file, ite_kernel_file, &
            ite_dt_adapt, ite_dt_max, ite_dt_grow, ite_dt_shrink)
       deallocate(rho_ite_init)
    end if
  end subroutine prepare_and_run_ite


  subroutine make_weights_from_existing_grid(n, r, w, dr_approx)
    integer, intent(in) :: n
    real(dp), intent(in) :: r(n)
    real(dp), intent(out) :: w(n), dr_approx
    integer :: i

    if (n < 3) stop 'ERROR: n must be >= 3 in make_weights_from_existing_grid.'
    w(1) = 0.5_dp * (r(2) - r(1))
    do i = 2, n-1
       w(i) = 0.5_dp * (r(i+1) - r(i-1))
    end do
    w(n) = 0.5_dp * (r(n) - r(n-1))

    ! The vW and ITE kinetic propagator currently assume a nearly uniform grid.
    dr_approx = (r(n) - r(1)) / real(n-1, dp)
  end subroutine make_weights_from_existing_grid

end program ks_skp_plus_nlf_constkf_v8
