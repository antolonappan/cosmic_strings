module cosmic_strings_cl_mod
  use iso_c_binding, only: c_double, c_int
  implicit none

  integer, parameter :: dp = c_double
  integer, parameter :: nz_bg = 1000

  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: c_tilde = 0.23_dp
  real(dp), parameter :: h0 = 70.0_dp
  real(dp), parameter :: h_hubble = h0 / 299792.458_dp
  real(dp), parameter :: omega_m = 0.3_dp
  real(dp), parameter :: omega_l = 0.7_dp
  real(dp), parameter :: chi_star = 14000.0_dp
  real(dp), parameter :: z_min = 0.01_dp
  real(dp), parameter :: z_max = 1100.0_dp

  logical, save :: bg_ready = .false.
  real(dp), save :: z_grid(nz_bg), chi_grid(nz_bg)

contains

  subroutine compute_cl_array_c(g_mu, p_rec, ell, n_ell, n_k, n_chi, k_min, k_max, cl) bind(C, name="compute_cl_array_c")
    real(c_double), value, intent(in) :: g_mu, p_rec, k_min, k_max
    integer(c_int), value, intent(in) :: n_ell, n_k, n_chi
    integer(c_int), intent(in) :: ell(n_ell)
    real(c_double), intent(out) :: cl(n_ell)

    integer :: i

    call init_background()

    do i = 1, n_ell
      cl(i) = compute_cl_one(int(ell(i)), g_mu, p_rec, int(n_k), int(n_chi), k_min, k_max)
    end do
  end subroutine compute_cl_array_c

  subroutine compute_cl_array(g_mu, p_rec, ell, n_k, n_chi, k_min, k_max, cl)
    real(dp), intent(in) :: g_mu, p_rec, k_min, k_max
    integer, intent(in) :: n_k, n_chi
    integer, intent(in) :: ell(:)
    real(dp), intent(out) :: cl(size(ell))

    integer :: i

    call init_background()

    do i = 1, size(ell)
      cl(i) = compute_cl_one(ell(i), g_mu, p_rec, n_k, n_chi, k_min, k_max)
    end do
  end subroutine compute_cl_array

  function compute_cl_one(ell, g_mu, p_rec, n_k, n_chi, k_min, k_max) result(cl)
    integer, intent(in) :: ell, n_k, n_chi
    real(dp), intent(in) :: g_mu, p_rec, k_min, k_max
    real(dp) :: cl

    real(dp), allocatable :: xk(:), wk(:), xchi(:), wchi(:)
    real(dp) :: log_k_min, log_k_max, log_k, k_val, inner, prefactor
    integer :: i

    allocate(xk(n_k), wk(n_k), xchi(n_chi), wchi(n_chi))

    log_k_min = log(k_min)
    log_k_max = log(k_max)
    call gauss_legendre(n_k, log_k_min, log_k_max, xk, wk)
    call gauss_legendre(n_chi, 0.0_dp, chi_star, xchi, wchi)

    cl = 0.0_dp
    do i = 1, n_k
      log_k = xk(i)
      k_val = exp(log_k)
      inner = inner_chi_integral(ell, k_val, g_mu, p_rec, xchi, wchi)
      cl = cl + wk(i) * inner * inner
    end do

    prefactor = 4.0_dp * pi / (real(ell, dp) * real(ell + 1, dp))
    cl = prefactor * cl

    deallocate(xk, wk, xchi, wchi)
  end function compute_cl_one

  function inner_chi_integral(ell, k_val, g_mu, p_rec, xchi, wchi) result(inner)
    integer, intent(in) :: ell
    real(dp), intent(in) :: k_val, g_mu, p_rec
    real(dp), intent(in) :: xchi(:), wchi(:)
    real(dp) :: inner

    real(dp) :: chi, x, d1, jl
    integer :: i

    inner = 0.0_dp
    do i = 1, size(xchi)
      chi = xchi(i)
      if (chi <= 0.0_dp) cycle
      x = k_val * chi
      jl = spherical_jn(ell, x)
      if (jl == 0.0_dp) cycle

      d1 = delta1sq(k_val, chi, g_mu, p_rec)
      if (d1 > 0.0_dp) then
        inner = inner + wchi(i) * sqrt(d1) * jl / chi
      end if
    end do
  end function inner_chi_integral

  function delta1sq(k_val, chi, g_mu, p_rec) result(delta)
    real(dp), intent(in) :: k_val, chi, g_mu, p_rec
    real(dp) :: delta

    real(dp) :: gamma, v_sq, v, xi, a, h_chi
    real(dp) :: prefactor, scale, er_term

    gamma = sqrt(pi * sqrt(2.0_dp) / (3.0_dp * c_tilde * p_rec))
    v_sq = 0.5_dp * (1.0_dp - pi / (3.0_dp * gamma))
    if (v_sq <= 0.0_dp) then
      delta = 0.0_dp
      return
    end if

    v = sqrt(v_sq)
    h_chi = h_of_chi(chi)
    a = a_of_chi(chi)
    xi = 1.0_dp / (h_chi * gamma)

    prefactor = (16.0_dp * pi * g_mu)**2
    prefactor = prefactor * (sqrt(6.0_dp * pi) * v**2 / (12.0_dp * (1.0_dp - v**2)))
    scale = (a / (k_val * xi))**5
    er_term = erf(k_val * xi / (a * 2.0_dp * sqrt(6.0_dp)))

    delta = prefactor * (4.0_dp * pi * k_val**3 * chi**2 * a**4 / h_chi) * scale * er_term
  end function delta1sq

  subroutine init_background()
    real(dp) :: dz, h_prev, h_cur
    integer :: i

    if (bg_ready) return

    dz = (z_max - z_min) / real(nz_bg - 1, dp)
    z_grid(1) = z_min
    chi_grid(1) = 0.0_dp
    h_prev = h_of_z(z_grid(1))

    do i = 2, nz_bg
      z_grid(i) = z_min + real(i - 1, dp) * dz
      h_cur = h_of_z(z_grid(i))
      chi_grid(i) = chi_grid(i - 1) + 0.5_dp * (1.0_dp / h_prev + 1.0_dp / h_cur) * dz
      h_prev = h_cur
    end do

    bg_ready = .true.
  end subroutine init_background

  function h_of_z(z) result(hz)
    real(dp), intent(in) :: z
    real(dp) :: hz

    hz = h_hubble * sqrt(omega_m * (1.0_dp + z)**3 + omega_l)
  end function h_of_z

  function z_of_chi(chi) result(z)
    real(dp), intent(in) :: chi
    real(dp) :: z

    integer :: lo, hi, mid
    real(dp) :: t

    call init_background()

    if (chi <= chi_grid(1)) then
      z = z_grid(1)
      return
    end if

    if (chi >= chi_grid(nz_bg)) then
      z = z_grid(nz_bg)
      return
    end if

    lo = 1
    hi = nz_bg
    do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if (chi_grid(mid) <= chi) then
        lo = mid
      else
        hi = mid
      end if
    end do

    t = (chi - chi_grid(lo)) / (chi_grid(hi) - chi_grid(lo))
    z = z_grid(lo) + t * (z_grid(hi) - z_grid(lo))
  end function z_of_chi

  function a_of_chi(chi) result(a)
    real(dp), intent(in) :: chi
    real(dp) :: a

    a = 1.0_dp / (1.0_dp + z_of_chi(chi))
  end function a_of_chi

  function h_of_chi(chi) result(hchi)
    real(dp), intent(in) :: chi
    real(dp) :: hchi

    hchi = h_of_z(z_of_chi(chi))
  end function h_of_chi

  function spherical_jn(ell, x) result(jl)
    integer, intent(in) :: ell
    real(dp), intent(in) :: x
    real(dp) :: jl

    real(dp) :: j0, j1, jm1, jcur, jp1
    real(dp) :: fnext, fcur, fprev, fl, f0, f1
    real(dp) :: true_j0, true_j1, norm, scale_guard
    integer :: n, m, extra, idx

    if (x == 0.0_dp) then
      if (ell == 0) then
        jl = 1.0_dp
      else
        jl = 0.0_dp
      end if
      return
    end if

    if (ell > 80 .and. x < 0.08_dp * real(ell, dp)) then
      jl = 0.0_dp
      return
    end if

    j0 = sin(x) / x
    if (ell == 0) then
      jl = j0
      return
    end if

    j1 = sin(x) / (x * x) - cos(x) / x
    if (ell == 1) then
      jl = j1
      return
    end if

    if (real(ell, dp) <= x) then
      jm1 = j0
      jcur = j1
      do n = 1, ell - 1
        jp1 = (real(2 * n + 1, dp) / x) * jcur - jm1
        jm1 = jcur
        jcur = jp1
      end do
      jl = jcur
      return
    end if

    extra = max(80, int(10.0_dp * sqrt(real(max(ell, 1), dp))))
    m = ell + extra

    fnext = 0.0_dp
    fcur = 1.0_dp
    fl = 0.0_dp
    f0 = 0.0_dp
    f1 = 0.0_dp

    do n = m, 1, -1
      fprev = (real(2 * n + 1, dp) / x) * fcur - fnext
      idx = n - 1

      if (idx == ell) fl = fprev
      if (idx == 1) f1 = fprev
      if (idx == 0) f0 = fprev

      fnext = fcur
      fcur = fprev

      scale_guard = max(abs(fnext), abs(fcur), abs(fl))
      if (scale_guard > 1.0e150_dp) then
        fnext = fnext * 1.0e-150_dp
        fcur = fcur * 1.0e-150_dp
        fl = fl * 1.0e-150_dp
        f0 = f0 * 1.0e-150_dp
        f1 = f1 * 1.0e-150_dp
      end if
    end do

    true_j0 = j0
    true_j1 = j1

    if (abs(true_j0) >= abs(true_j1) .and. abs(f0) > tiny(1.0_dp)) then
      norm = true_j0 / f0
    else if (abs(f1) > tiny(1.0_dp)) then
      norm = true_j1 / f1
    else
      norm = 0.0_dp
    end if

    jl = fl * norm
  end function spherical_jn

  subroutine gauss_legendre(n, a, b, x, w)
    integer, intent(in) :: n
    real(dp), intent(in) :: a, b
    real(dp), intent(out) :: x(n), w(n)

    integer :: i, j, m
    real(dp) :: xm, xl, z, z1, p1, p2, p3, pp
    real(dp), parameter :: eps = 1.0e-14_dp

    m = (n + 1) / 2
    xm = 0.5_dp * (b + a)
    xl = 0.5_dp * (b - a)

    do i = 1, m
      z = cos(pi * (real(i, dp) - 0.25_dp) / (real(n, dp) + 0.5_dp))

      do
        p1 = 1.0_dp
        p2 = 0.0_dp
        do j = 1, n
          p3 = p2
          p2 = p1
          p1 = ((2.0_dp * real(j, dp) - 1.0_dp) * z * p2 - (real(j, dp) - 1.0_dp) * p3) / real(j, dp)
        end do
        pp = real(n, dp) * (z * p1 - p2) / (z * z - 1.0_dp)
        z1 = z
        z = z1 - p1 / pp
        if (abs(z - z1) <= eps) exit
      end do

      x(i) = xm - xl * z
      x(n + 1 - i) = xm + xl * z
      w(i) = 2.0_dp * xl / ((1.0_dp - z * z) * pp * pp)
      w(n + 1 - i) = w(i)
    end do
  end subroutine gauss_legendre

end module cosmic_strings_cl_mod
