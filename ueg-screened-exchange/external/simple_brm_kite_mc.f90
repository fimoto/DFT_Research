program simple_brm_kite_mc
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 300)
  integer, parameter :: nquad = 64
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: alpha = (4.0_dp/(9.0_dp*pi))**(1.0_dp/3.0_dp)
  real(dp), parameter :: zeta3 = 1.20205690315959428539973816151145_dp
  real(dp), parameter :: onsager_ry = log(2.0_dp)/3.0_dp - 3.0_dp*zeta3/(2.0_dp*pi*pi)
  real(dp), parameter :: vball = 4.0_dp*pi/3.0_dp

  real(dp) :: rs, zeta, qscale
  integer :: nsample, nrep, qmap, proposal, seed0
  integer :: rep
  real(dp) :: gx(nquad), gw(nquad)
  real(dp) :: m1, m2, sum1, sum2, ss1, ss2
  real(dp) :: mu1, mu2, err1, err2, full, fullerr

  call parse_arguments(rs, zeta, nsample, nrep, qscale, qmap, proposal, seed0)
  call init_legendre_01(gx, gw)

  sum1 = 0.0_dp
  sum2 = 0.0_dp
  ss1 = 0.0_dp
  ss2 = 0.0_dp

  do rep = 1, nrep
     call set_seed(seed0 + 2003*rep + 11)
     m1 = term_mean(1, nsample, rs, zeta, qscale, qmap, proposal, gx, gw)
     call set_seed(seed0 + 2003*rep + 101)
     m2 = term_mean(2, nsample, rs, zeta, qscale, qmap, proposal, gx, gw)
     sum1 = sum1 + m1
     sum2 = sum2 + m2
     ss1 = ss1 + m1*m1
     ss2 = ss2 + m2*m2
  end do

  mu1 = sum1/real(nrep, dp)
  mu2 = sum2/real(nrep, dp)
  if (nrep > 1) then
     err1 = sqrt(max(0.0_dp, (ss1 - real(nrep,dp)*mu1*mu1)/real(nrep-1,dp)) &
          / real(nrep, dp))
     err2 = sqrt(max(0.0_dp, (ss2 - real(nrep,dp)*mu2*mu2)/real(nrep-1,dp)) &
          / real(nrep, dp))
  else
     err1 = 0.0_dp
     err2 = 0.0_dp
  end if
  full = onsager_ry + mu1 + mu2
  fullerr = sqrt(err1*err1 + err2*err2)

  write(*,'(a)') 'Simple BRM kite Monte Carlo'
  write(*,'(a,es14.6,a,f8.4)') '  rs = ', rs, ', zeta = ', zeta
  write(*,'(a,i0,a,i0)') '  samples/term/replicate = ', nsample, ', replicates = ', nrep
  write(*,'(a,f10.4,a,i0,a,i0)') '  q scale = ', qscale, ', q map = ', qmap, ', proposal = ', proposal
  write(*,'(a,es18.10,a,es10.2,a)') '  Delta1 = ', mu1, ' +- ', err1, ' Ry'
  write(*,'(a,es18.10,a,es10.2,a)') '  Delta2 = ', mu2, ' +- ', err2, ' Ry'
  write(*,'(a,f14.6,a,f10.6,a)') '  full kite = ', 1000.0_dp*full, ' +- ', 1000.0_dp*fullerr, ' mRy'

contains

  subroutine parse_arguments(rs, zeta, nsample, nrep, qscale, qmap, proposal, seed0)
    real(dp), intent(out) :: rs, zeta, qscale
    integer, intent(out) :: nsample, nrep, qmap, proposal, seed0
    character(len=128) :: arg
    integer :: narg

    rs = 5.0_dp
    zeta = 0.0_dp
    nsample = 4096
    nrep = 8
    qscale = 2.0_dp
    qmap = 1
    proposal = 1
    seed0 = 20260807

    narg = command_argument_count()
    if (narg >= 1) then
       call get_command_argument(1, arg); read(arg,*) rs
    end if
    if (narg >= 2) then
       call get_command_argument(2, arg); read(arg,*) zeta
    end if
    if (narg >= 3) then
       call get_command_argument(3, arg); read(arg,*) nsample
    end if
    if (narg >= 4) then
       call get_command_argument(4, arg); read(arg,*) nrep
    end if
    if (narg >= 5) then
       call get_command_argument(5, arg); read(arg,*) qscale
    end if
    if (narg >= 6) then
       call get_command_argument(6, arg); read(arg,*) qmap
    end if
    if (narg >= 7) then
       call get_command_argument(7, arg); read(arg,*) proposal
    end if
    if (narg >= 8) then
       call get_command_argument(8, arg); read(arg,*) seed0
    end if

    if (rs <= 0.0_dp) stop 'rs must be positive'
    if (abs(zeta) > 1.0_dp) stop 'zeta must be in [-1,1]'
    if (nsample < 1) stop 'nsample must be positive'
    if (nrep < 1) stop 'nrep must be positive'
    if (qscale <= 0.0_dp) stop 'qscale must be positive'
    if (qmap < 1 .or. qmap > 3) stop 'qmap must be 1 rational, 2 tangent, or 3 finite'
    if (proposal < 1 .or. proposal > 2) stop 'proposal must be 1 importance or 2 literal'
  end subroutine parse_arguments

  subroutine set_seed(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)
    call random_seed(size=n)
    allocate(put(n))
    do i = 1, n
       put(i) = modulo(seed + 104729*i + 7919*i*i, 2147483647)
       if (put(i) == 0) put(i) = i + 17
    end do
    call random_seed(put=put)
    deallocate(put)
  end subroutine set_seed

  real(dp) function term_mean(term, nsample, rs, zeta, qscale, qmap, proposal, gx, gw)
    integer, intent(in) :: term, nsample, qmap, proposal
    real(dp), intent(in) :: rs, zeta, qscale, gx(nquad), gw(nquad)
    integer :: i
    real(dp) :: u(9), val, sumv
    sumv = 0.0_dp
    do i = 1, nsample
       call random_number(u)
       if (proposal == 1) then
          val = sample_importance(term, u, rs, zeta, qscale, qmap, gx, gw)
       else
          val = sample_literal(term, u, rs, zeta, qscale, qmap, gx, gw)
       end if
       sumv = sumv + val
    end do
    term_mean = sumv/real(nsample, dp)
  end function term_mean

  real(dp) function sample_importance(term, uin, rs, zeta, qscale, qmap, gx, gw)
    integer, intent(in) :: term, qmap
    real(dp), intent(in) :: uin(9), rs, zeta, qscale, gx(nquad), gw(nquad)
    real(dp) :: u(9), lam, qv, qweight, k(3), q2(3), dir(3), center(3)
    real(dp) :: kvol, radius, rho, aq2, ak, tanu, pref, lfac, weight_common
    real(dp) :: fa, fb, dd, denom

    u = clip_unit(uin)
    lam = u(1)
    call map_q(u(2), qscale, qmap, qv, qweight)
    call sample_exit_cap(qv, u(3), u(4), u(5), k, kvol)
    call sample_direction(u(7), u(8), dir)

    if (term == 1) then
       center = k
       radius = 2.0_dp
    else
       center(1) = -k(1)
       center(2) = -k(2)
       center(3) = -k(3) - qv
       radius = sqrt(center(1)*center(1) + center(2)*center(2) + center(3)*center(3)) + 1.0_dp
    end if

    rho = radius*u(6)
    q2 = center + rho*dir
    if (.not. valid_exit_cap(q2, qv)) then
       sample_importance = 0.0_dp
       return
    end if

    aq2 = 0.5_dp*qv*qv + qv*q2(3)
    ak = 0.5_dp*qv*qv + qv*k(3)
    if (aq2 <= 0.0_dp .or. ak <= 0.0_dp) then
       sample_importance = 0.0_dp
       return
    end if

    tanu = tan(0.5_dp*pi*u(9))
    pref = 3.0_dp*alpha*rs/(16.0_dp*pi**7)
    lfac = pref*lam*lam/max(qv, 1.0e-300_dp)**4
    weight_common = qweight*kvol*(4.0_dp*pi*radius)*(0.5_dp*pi)*lfac

    if (term == 1) then
       dd = divided_difference(qv, aq2, ak, tanu, lam, rs, zeta, gx, gw)
       sample_importance = weight_common*dd
    else
       fa = f_spin_sum(qv, aq2*tanu, lam, rs, zeta, gx, gw)
       fb = f_spin_sum(qv, ak*tanu, lam, rs, zeta, gx, gw)
       denom = aq2 + ak
       sample_importance = -weight_common*(fa + fb)/denom
    end if
  end function sample_importance

  real(dp) function sample_literal(term, uin, rs, zeta, qscale, qmap, gx, gw)
    integer, intent(in) :: term, qmap
    real(dp), intent(in) :: uin(9), rs, zeta, qscale, gx(nquad), gw(nquad)
    real(dp) :: u(9), lam, qv, qweight, k(3), q2(3), kp(3), q2p(3), diff(3)
    real(dp) :: aq2, ak, tanu, pref, lfac, weight_common, fa, fb, dd, coul2

    u = clip_unit(uin)
    lam = u(1)
    call map_q(u(2), qscale, qmap, qv, qweight)
    call sample_unit_ball(u(3), u(4), u(5), k)
    call sample_unit_ball(u(6), u(7), u(8), q2)

    kp = k
    q2p = q2
    kp(3) = kp(3) + qv
    q2p(3) = q2p(3) + qv
    if (dot3(kp,kp) <= 1.0_dp .or. dot3(q2p,q2p) <= 1.0_dp) then
       sample_literal = 0.0_dp
       return
    end if

    aq2 = 0.5_dp*qv*qv + qv*q2(3)
    ak = 0.5_dp*qv*qv + qv*k(3)
    if (aq2 <= 0.0_dp .or. ak <= 0.0_dp) then
       sample_literal = 0.0_dp
       return
    end if

    tanu = tan(0.5_dp*pi*u(9))
    pref = 3.0_dp*alpha*rs/(16.0_dp*pi**7)
    lfac = pref*lam*lam/max(qv, 1.0e-300_dp)**4
    weight_common = qweight*vball*vball*(0.5_dp*pi)*lfac

    if (term == 1) then
       diff = k - q2
       coul2 = max(dot3(diff,diff), 1.0e-300_dp)
       dd = divided_difference(qv, aq2, ak, tanu, lam, rs, zeta, gx, gw)
       sample_literal = weight_common*dd/coul2
    else
       diff(1) = k(1) + q2(1)
       diff(2) = k(2) + q2(2)
       diff(3) = k(3) + q2(3) + qv
       coul2 = max(dot3(diff,diff), 1.0e-300_dp)
       fa = f_spin_sum(qv, aq2*tanu, lam, rs, zeta, gx, gw)
       fb = f_spin_sum(qv, ak*tanu, lam, rs, zeta, gx, gw)
       sample_literal = -weight_common*(fa + fb)/(coul2*(aq2 + ak))
    end if
  end function sample_literal

  subroutine map_q(tin, scale, qmap, qv, qweight)
    real(dp), intent(in) :: tin, scale
    integer, intent(in) :: qmap
    real(dp), intent(out) :: qv, qweight
    real(dp) :: t, angle, dqdt
    t = max(1.0e-12_dp, min(1.0_dp - 1.0e-12_dp, tin))
    if (qmap == 1) then
       qv = scale*t/(1.0_dp - t)
       dqdt = scale/(1.0_dp - t)**2
    else if (qmap == 2) then
       angle = 0.5_dp*pi*t
       qv = scale*tan(angle)
       dqdt = 0.5_dp*pi*scale/(cos(angle)**2)
    else
       qv = scale*t
       dqdt = scale
    end if
    qweight = 4.0_dp*pi*qv*qv*dqdt
  end subroutine map_q

  subroutine sample_unit_ball(u1, u2, u3, p)
    real(dp), intent(in) :: u1, u2, u3
    real(dp), intent(out) :: p(3)
    real(dp) :: r, ct, st, ph
    r = max(1.0e-12_dp, min(1.0_dp - 1.0e-12_dp, u1))**(1.0_dp/3.0_dp)
    ct = 2.0_dp*u2 - 1.0_dp
    ph = 2.0_dp*pi*u3
    st = sqrt(max(0.0_dp, 1.0_dp - ct*ct))
    p(1) = r*st*cos(ph)
    p(2) = r*st*sin(ph)
    p(3) = r*ct
  end subroutine sample_unit_ball

  subroutine sample_direction(u1, u2, d)
    real(dp), intent(in) :: u1, u2
    real(dp), intent(out) :: d(3)
    real(dp) :: ct, st, ph
    ct = 2.0_dp*u1 - 1.0_dp
    ph = 2.0_dp*pi*u2
    st = sqrt(max(0.0_dp, 1.0_dp - ct*ct))
    d(1) = st*cos(ph)
    d(2) = st*sin(ph)
    d(3) = ct
  end subroutine sample_direction

  subroutine sample_exit_cap(q, u1, u2, u3, p, vol)
    real(dp), intent(in) :: q, u1, u2, u3
    real(dp), intent(out) :: p(3), vol
    real(dp) :: prob1, v1, z, rho2, rho, ph, t, z0, f0, f1, target
    real(dp) :: lo, hi, mid, fm, upper2, lower2
    integer :: it

    vol = exit_cap_volume(q)
    if (q >= 2.0_dp) then
       call sample_unit_ball(u1, u2, u3, p)
       return
    end if

    v1 = pi*q*(1.0_dp - 0.5_dp*q)**2
    prob1 = max(0.0_dp, min(1.0_dp, v1/max(vol, 1.0e-300_dp)))
    if (u1 < prob1) then
       t = u1/max(prob1, 1.0e-300_dp)
       z = -0.5_dp*q + (1.0_dp - 0.5_dp*q)*sqrt(t)
    else
       z0 = 1.0_dp - q
       t = (u1 - prob1)/max(1.0_dp - prob1, 1.0e-300_dp)
       f0 = z0 - z0**3/3.0_dp
       f1 = 2.0_dp/3.0_dp
       target = f0 + t*(f1 - f0)
       lo = z0
       hi = 1.0_dp
       do it = 1, 50
          mid = 0.5_dp*(lo + hi)
          fm = mid - mid**3/3.0_dp
          if (fm < target) then
             lo = mid
          else
             hi = mid
          end if
       end do
       z = 0.5_dp*(lo + hi)
    end if
    upper2 = max(0.0_dp, 1.0_dp - z*z)
    lower2 = max(0.0_dp, 1.0_dp - (z + q)**2)
    rho2 = lower2 + u2*max(0.0_dp, upper2 - lower2)
    rho = sqrt(rho2)
    ph = 2.0_dp*pi*u3
    p(1) = rho*cos(ph)
    p(2) = rho*sin(ph)
    p(3) = z
  end subroutine sample_exit_cap

  real(dp) function exit_cap_volume(q)
    real(dp), intent(in) :: q
    if (q >= 2.0_dp) then
       exit_cap_volume = vball
    else
       exit_cap_volume = pi*q*(12.0_dp - q*q)/12.0_dp
    end if
  end function exit_cap_volume

  logical function valid_exit_cap(p, q)
    real(dp), intent(in) :: p(3), q
    real(dp) :: ps(3)
    ps = p
    ps(3) = ps(3) + q
    valid_exit_cap = (dot3(p,p) <= 1.0_dp) .and. (dot3(ps,ps) > 1.0_dp)
  end function valid_exit_cap

  real(dp) function divided_difference(q, a, b, tanu, lam, rs, zeta, gx, gw)
    real(dp), intent(in) :: q, a, b, tanu, lam, rs, zeta, gx(nquad), gw(nquad)
    real(dp) :: fa, fb, den, scale, nm
    fa = f_spin_sum(q, a*tanu, lam, rs, zeta, gx, gw)
    fb = f_spin_sum(q, b*tanu, lam, rs, zeta, gx, gw)
    den = a - b
    scale = max(1.0_dp, abs(a) + abs(b))
    if (abs(den) < 2.0e-7_dp*scale) then
       nm = 0.5_dp*(a*tanu + b*tanu)
       divided_difference = tanu*df_spin_sum(q, nm, lam, rs, zeta, gx, gw)
    else
       divided_difference = (fa - fb)/den
    end if
  end function divided_difference

  real(dp) function f_spin_sum(q, nu, lam, rs, zeta, gx, gw)
    real(dp), intent(in) :: q, nu, lam, rs, zeta, gx(nquad), gw(nquad)
    f_spin_sum = f_lambda_sigma(q, nu, lam, rs, zeta, 1, gx, gw) &
               + f_lambda_sigma(q, nu, lam, rs, zeta, -1, gx, gw)
  end function f_spin_sum

  real(dp) function f_lambda_sigma(q, nu, lam, rs, zeta, sigma, gx, gw)
    real(dp), intent(in) :: q, nu, lam, rs, zeta, gx(nquad), gw(nquad)
    integer, intent(in) :: sigma
    real(dp) :: xu, xd, x, xp, theta, ssum, den

    call spin_radii(zeta, xu, xd)
    if (sigma > 0) then
       x = xu
    else
       x = xd
    end if
    if (x < 1.0e-14_dp) then
       f_lambda_sigma = 0.0_dp
       return
    end if

    ssum = 0.0_dp
    xp = xu
    if (xp > 1.0e-14_dp) then
       theta = x/xp
       ssum = ssum + xp*g_brm_quad(theta*q, theta*theta*abs(nu), gx, gw)
    end if
    xp = xd
    if (xp > 1.0e-14_dp) then
       theta = x/xp
       ssum = ssum + xp*g_brm_quad(theta*q, theta*theta*abs(nu), gx, gw)
    end if

    den = x*x + lam*alpha*rs*ssum/(pi*max(q, 1.0e-300_dp)**2)
    f_lambda_sigma = x**3*ssum/den
  end function f_lambda_sigma

  real(dp) function df_spin_sum(q, nu, lam, rs, zeta, gx, gw)
    real(dp), intent(in) :: q, nu, lam, rs, zeta, gx(nquad), gw(nquad)
    real(dp) :: h, lo, hi
    if (abs(nu) < 1.0e-12_dp) then
       df_spin_sum = 0.0_dp
       return
    end if
    h = 2.0e-5_dp*max(1.0_dp, abs(nu))
    lo = max(0.0_dp, abs(nu) - h)
    hi = abs(nu) + h
    df_spin_sum = (f_spin_sum(q, hi, lam, rs, zeta, gx, gw) &
                 - f_spin_sum(q, lo, lam, rs, zeta, gx, gw))/(hi - lo)
  end function df_spin_sum

  real(dp) function g_brm_quad(q, nu, gx, gw)
    real(dp), intent(in) :: q, nu, gx(nquad), gw(nquad)
    integer :: i
    real(dp) :: kk, a, den, inc, qq, nn
    qq = max(abs(q), 1.0e-12_dp)
    nn = abs(nu)
    a = 0.5_dp*qq*qq
    g_brm_quad = 0.0_dp
    do i = 1, nquad
       kk = gx(i)
       den = (qq*kk - a)**2 + nn*nn
       den = max(den, 1.0e-300_dp)
       inc = log(1.0_dp + 4.0_dp*a*qq*kk/den)
       g_brm_quad = g_brm_quad + gw(i)*kk*inc
    end do
    g_brm_quad = g_brm_quad/qq
  end function g_brm_quad

  subroutine spin_radii(zeta, xu, xd)
    real(dp), intent(in) :: zeta
    real(dp), intent(out) :: xu, xd
    xu = max(0.0_dp, 1.0_dp + zeta)**(1.0_dp/3.0_dp)
    xd = max(0.0_dp, 1.0_dp - zeta)**(1.0_dp/3.0_dp)
  end subroutine spin_radii

  real(dp) function dot3(a, b)
    real(dp), intent(in) :: a(3), b(3)
    dot3 = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
  end function dot3

  elemental real(dp) function clip_scalar(x)
    real(dp), intent(in) :: x
    clip_scalar = max(1.0e-12_dp, min(1.0_dp - 1.0e-12_dp, x))
  end function clip_scalar

  elemental function clip_unit(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y = clip_scalar(x)
  end function clip_unit

  subroutine init_legendre_01(x, w)
    real(dp), intent(out) :: x(nquad), w(nquad)
    integer :: i, j, m
    real(dp) :: z, z1, p1, p2, p3, pp
    m = (nquad + 1)/2
    do i = 1, m
       z = cos(pi*(real(i,dp) - 0.25_dp)/(real(nquad,dp) + 0.5_dp))
       do
          p1 = 1.0_dp
          p2 = 0.0_dp
          do j = 1, nquad
             p3 = p2
             p2 = p1
             p1 = ((2.0_dp*real(j,dp) - 1.0_dp)*z*p2 - (real(j,dp) - 1.0_dp)*p3)/real(j,dp)
          end do
          pp = real(nquad,dp)*(z*p1 - p2)/(z*z - 1.0_dp)
          z1 = z
          z = z1 - p1/pp
          if (abs(z - z1) < 1.0e-15_dp) exit
       end do
       x(i) = 0.5_dp*(1.0_dp - z)
       x(nquad + 1 - i) = 0.5_dp*(1.0_dp + z)
       w(i) = 1.0_dp/((1.0_dp - z*z)*pp*pp)
       w(nquad + 1 - i) = w(i)
    end do
  end subroutine init_legendre_01

end program simple_brm_kite_mc
