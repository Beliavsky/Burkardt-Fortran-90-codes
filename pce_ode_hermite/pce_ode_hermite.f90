function he_double_product_integral ( i, j )

!*****************************************************************************80
!
!! he_double_product_integral(): integral of He(i,x)*He(j,x)*e^(-x^2/2).
!
!  Discussion:
!
!    VALUE = integral ( -oo < x < +oo ) He(i,x)*He(j,x) exp(-x^2/2) dx
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical Methods for Stochastic Computations: A Spectral Method Approach,
!    Princeton, 2010,
!    ISBN13: 978-0-691-14212-8,
!    LC: QA274.23.X58.
!
!  Parameters:
!
!    Input, integer I, J, the polynomial indices.
!
!    Output, real ( kind = rk ) HE_DOUBLE_PRODUCT_INTEGRAL, the value of 
!    the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) he_double_product_integral
  integer i
  integer j
  real ( kind = rk ) r8_factorial
  real ( kind = rk ) value

  if ( i == j ) then
    value = r8_factorial ( i )
  else
    value = 0.0D+00
  end if

  he_double_product_integral = value

  return
end
function he_triple_product_integral ( i, j, k )

!*****************************************************************************80
!
!! he_triple_product_integral: integral of He(i,x)*He(j,x)*He(k,x)*e^(-x^2/2).
!
!  Discussion:
!
!    VALUE = integral ( -oo < x < +oo ) He(i,x)*He(j,x)*He(k,x) exp(-x^2/2) dx
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical Methods for Stochastic Computations: A Spectral Method Approach,
!    Princeton, 2010,
!    ISBN13: 978-0-691-14212-8,
!    LC: QA274.23.X58.
!
!  Parameters:
!
!    Input, integer I, J, K, the polynomial indices.
!
!    Output, real ( kind = rk ) HE_TRIPLE_PRODUCT_INTEGRAL, the value of the
!    integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) he_triple_product_integral
  integer i
  integer j
  integer k
  real ( kind = rk ) r8_factorial
  integer s
  real ( kind = rk ) value

  s = ( i + j + k ) / 2

  if ( s < max ( i, j, k ) ) then
    value = 0.0D+00
  else if ( mod ( i + j + k, 2 ) /= 0 ) then
    value = 0.0D+00
  else
    value = r8_factorial ( i ) / r8_factorial ( s - i ) &
          * r8_factorial ( j ) / r8_factorial ( s - j ) &
          * r8_factorial ( k ) / r8_factorial ( s - k )
  end if

  he_triple_product_integral = value

  return
end
subroutine pce_ode_hermite ( ti, tf, nt, ui, np, alpha_mu, alpha_sigma, t, u )

!*****************************************************************************80
!
!! pce_ode_hermite() applies the polynomial chaos expansion to a scalar ODE.
!
!  Discussion:
!
!    The deterministic equation is
!
!      du/dt = - alpha * u,
!      u(0) = u0
!
!    In the stochastic version, it is assumed that the decay coefficient
!    ALPHA is a Gaussian random variable with mean value ALPHA_MU and variance
!    ALPHA_SIGMA^2.
!
!    The exact expected value of the stochastic equation will be
!
!      u(t) = u0 * exp ( t^2/2)
!
!    This should be matched by the first component of the polynomial chaos
!    expansion.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) TI, TF, the initial and final times.
!
!    Input, integer NT, the number of output points.
!
!    Input, real ( kind = rk ) UI, the initial condition.
!
!    Input, integer NP, the degree of the expansion.  Polynomials 
!    of degree 0 through NP will be used.
!
!    Input, real ( kind = rk ) ALPHA_MU, ALPHA_SIGMA, the mean and standard 
!    deviation of the decay coefficient.
!
!    Output, real ( kind = rk ) T(0:NT), U(0:NT,0:NP), the times and the PCE 
!    coefficients at the successive time steps.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer np
  integer nt

  real ( kind = rk ) alpha_mu
  real ( kind = rk ) alpha_sigma
  real ( kind = rk ) dp
  real ( kind = rk ) dt
  real ( kind = rk ) he_double_product_integral
  real ( kind = rk ) he_triple_product_integral
  integer i
  integer it
  integer j
  integer k
  real ( kind = rk ) t(0:nt)
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) term
  real ( kind = rk ) tf
  real ( kind = rk ) ti
  real ( kind = rk ) tp
  real ( kind = rk ) u(0:nt,0:np)
  real ( kind = rk ) u1(0:np)
  real ( kind = rk ) u2(0:np)
  real ( kind = rk ) ui

  dt = ( tf - ti ) / real ( nt, kind = rk )
!
!  Set the PCE coefficients for the initial time.
!
  t1 = ti
  u1(0) = ui
  u1(1:np) = 0.0D+00
!
!  Copy into the output arrays.
!
  t(0) = t1
  u(0,0:np) = u1(0:np)
!
!  Time integration.
!
  do it = 1, nt

    t2 = ( real ( nt - it, kind = rk ) * ti   &
         + real (      it, kind = rk ) * tf ) &
         / real ( nt,      kind = rk )

    do k = 0, np

      dp = he_double_product_integral ( k, k )

      term = - alpha_mu * u1(k)

      i = 1
      do j = 0, np
        tp = he_triple_product_integral ( i, j, k )
        term = term - alpha_sigma * u1(j) * tp / dp
      end do

      u2(k) = u1(k) + dt * term

    end do
!
!  Prepare for next step.
!
    t1 = t2
    u1(0:np) = u2(0:np)
!
!  Copy into the output arrays.
!
    t(it) = t1
    u(it,0:np) = u1(0:np)

  end do

  return
end
function r8_factorial ( n )

!*****************************************************************************80
!
!! r8_factorial() computes the factorial of N.
!
!  Discussion:
!
!    factorial ( N ) = product ( 1 <= I <= N ) I
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the argument of the factorial function.
!    If N is less than 1, the function value is returned as 1.
!
!    Output, real ( kind = rk ) R8_FACTORIAL, the factorial of N.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_factorial
  integer i
  integer n

  r8_factorial = 1.0D+00

  do i = 1, n
    r8_factorial = r8_factorial * real ( i, kind = rk )
  end do

  return
end

