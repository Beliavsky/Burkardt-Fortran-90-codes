subroutine gegenbauer_alpha_check ( alpha, check )

!*****************************************************************************80
!
!! gegenbauer_alpha_check() checks the value of ALPHA.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 November 2015
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) ALPHA, a parameter which is part of the 
!    definition of the Gegenbauer polynomials.  It must be greater than -0.5.
!
!  Output:
!
!    logical CHECK.
!    TRUE, ALPHA is acceptable.
!    FALSE, ALPHA is not acceptable. 
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alpha
  logical check
  logical squawk

  squawk = .false.

  if ( - 0.5D+00 < alpha ) then


    check = .true.

  else

    check = .false.

    if ( squawk ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'gegenbauer_polynomial_value(): Fatal error!';
      write ( *, '(a)' ) '  Illegal value of ALPHA.'
      write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha 
      write ( *, '(a)' ) '  but ALPHA must be greater than -0.5.'
    end if

  end if

  return
end
subroutine gegenbauer_ek_compute ( n, alpha, x, w )

!*****************************************************************************80
!
!! gegenbauer_ek_compute() computes a Gauss-Gegenbauer quadrature rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= 1 ) (1-x^2)^(alpha-1/2) * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 November 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Input:
!
!    integer N, the order.
!
!    real ( kind = rk ) ALPHA, the exponent of (1-X*X) in the weight.
!    -1.0 < ALPHA.
!
!  Output:
!
!    real ( kind = rk ) X(N), the abscissas.
!
!    real ( kind = rk ) W(N), the weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) abi
  real ( kind = rk ) alpha
  real ( kind = rk ) bj(n)
  integer i
  real ( kind = rk ) i_r8
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) zemu
!
!  Check N.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'gegenbauer_ek_compute(): Fatal error!'
    write ( *, '(a)' ) '  1 <= N is required.'
    stop 1
  end if
!
!  Check ALPHA.
!
  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'gegenbauer_ek_compute(): Fatal error!'
    write ( *, '(a)' ) '  -1.0 < ALPHA is required.'
    stop 1
  end if
!
!  Define the zero-th moment.
!
  zemu = 2.0 ** ( 2.0 * alpha + 1.0 ) &
    * gamma ( alpha + 1.0 ) &
    * gamma ( alpha + 1.0 ) &
    / gamma ( 2.0 * alpha + 2.0 )
!
!  Define the Jacobi matrix.
!
  x(1:n) = 0.0D+00

  bj(1) = 4.0D+00 * ( alpha + 1.0D+00 ) ** 2 &
    / ( ( 2.0D+00 * alpha + 3.0D+00 ) * ( 2.0D+00 * alpha + 2.0D+00 ) ** 2 )

  do i = 2, n
    i_r8 = real ( i, kind = rk )
    abi = 2.0D+00 * ( alpha + i_r8 );
    bj(i) = 4.0D+00 * i_r8 * ( alpha + i_r8 ) ** 2 &
      * ( 2.0D+00 * alpha + i_r8 ) &
      / ( ( abi - 1.0D+00 ) * ( abi + 1.0D+00 ) * abi * abi )
  end do

  bj(1:n) = sqrt ( bj(1:n) )

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n) ** 2

  return
end
subroutine gegenbauer_integral ( expon, alpha, value )

!*****************************************************************************80
!
!! gegenbauer_integral(): integral of a monomial with Gegenbauer weight.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= +1 ) x^expon (1-x*x)^(alpha-1/2) dx
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer EXPON, the exponent.
!
!    real ( kind = rk ) ALPHA, part of the exponent of (1-X^2).
!    -0.5 < ALPHA.
!
!  Output:
!
!    real ( kind = rk ) VALUE, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alpha
  real ( kind = rk ) arg1
  real ( kind = rk ) arg2
  real ( kind = rk ) arg3
  real ( kind = rk ) arg4
  real ( kind = rk ) c
  integer expon
  real ( kind = rk ) value
  real ( kind = rk ) value1

  if ( mod ( expon, 2 ) == 1 ) then
    value = 0.0D+00
    return
  end if

  c = real ( expon, kind = rk )

  arg1 = - alpha
  arg2 =   1.0D+00 + c
  arg3 =   2.0D+00 + alpha + c
  arg4 = - 1.0D+00

  call r8_hyper_2f1 ( arg1, arg2, arg3, arg4, value1 )

  value = 2.0D+00 * gamma ( 1.0D+00 + c ) * gamma ( 1.0D+00 + alpha ) &
    * value1 / gamma ( 2.0D+00 + alpha + c )

  return
end
subroutine gegenbauer_polynomial_value ( m, n, alpha, x, c )

!*****************************************************************************80
!
!! gegenbauer_polynomial_value() computes the Gegenbauer polynomials C(I,ALPHA,X).
!
!  Discussion:
!
!    The Gegenbauer polynomial can be evaluated in Mathematica with
!    the command 
!
!      GegenbauerC[n,m,x]
!
!    ALPHA must be greater than -0.5.
!
!    If ALPHA = 1, the Gegenbauer polynomials reduce to the Chebyshev
!    polynomials of the second kind.
!
!  Differential equation:
!
!    (1-X*X) Y'' - (2 ALPHA + 1) X Y' + N (N + 2 ALPHA) Y = 0
!
!  Recursion:
!
!    C(0,ALPHA,X) = 1,
!    C(1,ALPHA,X) = 2*ALPHA*X
!    C(N,ALPHA,X) = ( (2*N-2+2*ALPHA) * X * C(N-1,ALPHA,X) 
!                   + ( -N+2-2*ALPHA)     * C(N-2,ALPHA,X) ) / N
!
!  Norm:
!
!    Integral ( -1 <= X <= 1 ) 
!      ( 1 - X^2 )^( ALPHA - 0.5 ) * C(N,ALPHA,X)^2 dX
!
!    = PI * 2^( 1 - 2 * ALPHA ) * Gamma ( N + 2 * ALPHA ) 
!      / ( N! * ( N + ALPHA ) * ( Gamma ( ALPHA ) )^2 )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 November 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Input:
!
!    integer M, the highest order polynomial to compute.
!    Note that polynomials 0 through M will be computed.
!
!    integer N, the number of evaluation points.
!
!    real ( kind = rk ) ALPHA, a parameter which is part of the 
!    definition of the Gegenbauer polynomials.  It must be greater than -0.5.
!
!    real ( kind = rk ) X(N), the evalulation points.
!
!  Output:
!
!    real ( kind = rk ) C(0:M,N), the values of the first N+1 Gegenbauer
!    polynomials at the point X.  
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) alpha
  real ( kind = rk ) c(0:m,n)
  integer i
  real ( kind = rk ) x(n)

  if ( alpha <= -0.5D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'gegenbauer_polynomial_value(): Fatal error!'
    write ( *, '(a,g14.6)' ) '  Illegal value of ALPHA = ', alpha
    write ( *, '(a)' ) '  but ALPHA must be greater than -0.5.'
    stop 1
  end if

  if ( m < 0 ) then
    return
  end if

  if ( n < 1 ) then
    return
  end if

  c(0,1:n) = 1.0D+00

  if ( m < 1 ) then
    return
  end if

  c(1,1:n) = 2.0D+00 * alpha * x(1:n)

  do i = 2, m
    c(i,1:n) = &
      ( ( real ( 2 * i - 2, kind = rk ) + 2.0D+00 * alpha ) * x(1:n) * c(i-1,1:n)   &
      + ( real (   - i + 2, kind = rk ) - 2.0D+00 * alpha )          * c(i-2,1:n) ) &
      /   real (     i,     kind = rk )
  end do
 
  return
end
subroutine gegenbauer_polynomial_values ( n_data, n, a, x, fx )

!*****************************************************************************80
!
!! gegenbauer_polynomial_values() returns some values of Gegenbauer polynomials.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Input:
!
!    integer N_DATA.  The user sets N_DATA to 0 before the first call.  
!
!  Output:
!
!    integer N_DATA.  The routine increments N_DATA by 1, 
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    integer N, the order parameter of the function.
!
!    real ( kind = rk ) A, the real parameter of the function.
!
!    real ( kind = rk ) X, the argument of the function.
!
!    real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 38

  real ( kind = rk ) a
  real ( kind = rk ), save, dimension ( n_max ) :: a_vec = (/ &
     0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.0D+00, &
     1.0D+00,  2.0D+00,  3.0D+00, &
     4.0D+00,  5.0D+00,  6.0D+00, &
     7.0D+00,  8.0D+00,  9.0D+00, &
    10.0D+00,  3.0D+00,  3.0D+00, &
     3.0D+00,  3.0D+00,  3.0D+00, &
     3.0D+00,  3.0D+00,  3.0D+00, &
     3.0D+00,  3.0D+00,  3.0D+00, &
     3.0D+00,  3.0D+00,  3.0D+00, &
     3.0D+00,  3.0D+00 /)
  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
    1.0000000000D+00,   0.2000000000D+00,  -0.4400000000D+00, &
   -0.2800000000D+00,   0.2320000000D+00,   0.3075200000D+00, &
   -0.0805760000D+00,  -0.2935168000D+00,  -0.0395648000D+00, &
    0.2459712000D+00,   0.1290720256D+00,   0.0000000000D+00, &
   -0.3600000000D+00,  -0.0800000000D+00,   0.8400000000D+00, &
    2.4000000000D+00,   4.6000000000D+00,   7.4400000000D+00, &
   10.9200000000D+00,  15.0400000000D+00,  19.8000000000D+00, &
   25.2000000000D+00,  -9.0000000000D+00,  -0.1612800000D+00, &
   -6.6729600000D+00,  -8.3750400000D+00,  -5.5267200000D+00, &
    0.0000000000D+00,   5.5267200000D+00,   8.3750400000D+00, &
    6.6729600000D+00,   0.1612800000D+00,  -9.0000000000D+00, &
  -15.4252800000D+00,  -9.6969600000D+00,  22.4409600000D+00, &
  100.8892800000D+00, 252.0000000000D+00 /)
  integer n
  integer n_data
  integer, save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10,  2, &
     2,  2,  2, &
     2,  2,  2, &
     2,  2,  2, &
     2,  5,  5, &
     5,  5,  5, &
     5,  5,  5, &
     5,  5,  5, &
     5,  5,  5, &
     5,  5 /)
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( n_max ) :: x_vec = (/ &
    0.20D+00,  0.20D+00,  0.20D+00, &
    0.20D+00,  0.20D+00,  0.20D+00, &
    0.20D+00,  0.20D+00,  0.20D+00, &
    0.20D+00,  0.20D+00,  0.40D+00, &
    0.40D+00,  0.40D+00,  0.40D+00, &
    0.40D+00,  0.40D+00,  0.40D+00, &
    0.40D+00,  0.40D+00,  0.40D+00, &
    0.40D+00, -0.50D+00, -0.40D+00, &
   -0.30D+00, -0.20D+00, -0.10D+00, &
    0.00D+00,  0.10D+00,  0.20D+00, &
    0.30D+00,  0.40D+00,  0.50D+00, &
    0.60D+00,  0.70D+00,  0.80D+00, &
    0.90D+00,  1.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    a = 0.0
    x = 0.0
    fx = 0.0
  else
    n = n_vec(n_data)
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine gegenbauer_ss_compute ( n, alpha, x, w )

!*****************************************************************************80
!
!! gegenbauer_ss_compute() computes a Gauss-Gegenbauer quadrature rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= 1 ) (1-x^2)^(alpha-1/2) * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    Thanks to Janiki Raman for pointing out a problem in an earlier
!    version of the code that occurred when ALPHA was -0.5.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Input:
!
!    integer N, the order.
!
!    real ( kind = rk ) ALPHA, the exponent of (1-X*X) in the weight.
!    -1.0 < ALPHA.
!
!  Output:
!
!    real ( kind = rk ) X(N), the abscissas.
!
!    real ( kind = rk ) W(N), the weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) alpha
  real ( kind = rk ) an
  real ( kind = rk ) c(n)
  real ( kind = rk ) cc
  real ( kind = rk ) delta
  real ( kind = rk ) dp2
  integer i
  real ( kind = rk ) p1
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) r3
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) xval
!
!  Check N.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'gegenbauer_ss_compute(): Fatal error!'
    write ( *, '(a)' ) '  1 <= N is required.'
    stop 1
  end if
!
!  Check ALPHA.
!
  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'gegenbauer_ss_compute(): Fatal error!'
    write ( *, '(a)' ) '  -1.0 < ALPHA is required.'
    stop 1
  end if
!
!  Set the recursion coefficients.
!
  c(1) = 0.0D+00

  if ( 2 <= n ) then
    c(2) = 1.0D+00 / ( 2.0D+00 * alpha + 3.0D+00 )
  end if

  do i = 3, n

    c(i) = real ( i - 1, kind = rk ) &
          * ( alpha + alpha + real ( i - 1, kind = rk ) ) / &
          ( ( alpha + alpha + real ( 2 * i - 1, kind = rk ) ) &
          * ( alpha + alpha + real ( 2 * i - 3, kind = rk ) ) )

  end do

  delta = gamma ( alpha         + 1.0D+00 ) &
        * gamma (         alpha + 1.0D+00 ) &
        / gamma ( alpha + alpha + 2.0D+00 )

  cc = delta * 2.0D+00 ** ( 2.0D+00 * alpha + 1.0D+00 ) * product ( c(2:n) )

  do i = 1, n

    if ( i == 1 ) then

      an = alpha / real ( n, kind = rk )

      r1 = ( 1.0D+00 + alpha ) &
        * ( 2.78D+00 / ( 4.0D+00 + real ( n * n, kind = rk ) ) &
        + 0.768D+00 * an / real ( n, kind = rk ) )

      r2 = 1.0D+00 + 2.44D+00 * an + 1.282D+00 * an * an

      xval = ( r2 - r1 ) / r2

    else if ( i == 2 ) then

      r1 = ( 4.1D+00 + alpha ) / &
        ( ( 1.0D+00 + alpha ) * ( 1.0D+00 + 0.156D+00 * alpha ) )

      r2 = 1.0D+00 + 0.06D+00 * ( real ( n, kind = rk ) - 8.0D+00 ) * &
        ( 1.0D+00 + 0.12D+00 * alpha ) / real ( n, kind = rk )

      r3 = 1.0D+00 + 0.012D+00 * alpha * &
        ( 1.0D+00 + 0.25D+00 * abs ( alpha ) ) / real ( n, kind = rk )

      xval = xval - r1 * r2 * r3 * ( 1.0D+00 - xval )

    else if ( i == 3 ) then

      r1 = ( 1.67D+00 + 0.28D+00 * alpha ) / ( 1.0D+00 + 0.37D+00 * alpha )

      r2 = 1.0D+00 + 0.22D+00 * ( real ( n, kind = rk ) - 8.0D+00 ) &
        / real ( n, kind = rk )

      r3 = 1.0D+00 + 8.0D+00 * alpha / &
        ( ( 6.28D+00 + alpha ) * real ( n * n, kind = rk ) )

      xval = xval - r1 * r2 * r3 * ( x(1) - xval )

    else if ( i < n - 1 ) then

      xval = 3.0D+00 * x(i-1) - 3.0D+00 * x(i-2) + x(i-3)

    else if ( i == n - 1 ) then

      r1 = ( 1.0D+00 + 0.235D+00 * alpha ) / ( 0.766D+00 + 0.119D+00 * alpha )

      r2 = 1.0D+00 / ( 1.0D+00 + 0.639D+00 &
        * ( real ( n, kind = rk ) - 4.0D+00 ) &
        / ( 1.0D+00 + 0.71D+00 * ( real ( n, kind = rk ) - 4.0D+00 ) ) )

      r3 = 1.0D+00 / ( 1.0D+00 + 20.0D+00 * alpha / ( ( 7.5D+00 + alpha ) * &
        real ( n * n, kind = rk ) ) )

      xval = xval + r1 * r2 * r3 * ( xval - x(i-2) )

    else if ( i == n ) then

      r1 = ( 1.0D+00 + 0.37D+00 * alpha ) / ( 1.67D+00 + 0.28D+00 * alpha )

      r2 = 1.0D+00 / &
        ( 1.0D+00 + 0.22D+00 * ( real ( n, kind = rk ) - 8.0D+00 ) &
        / real ( n, kind = rk ) )

      r3 = 1.0D+00 / ( 1.0D+00 + 8.0D+00 * alpha / &
        ( ( 6.28D+00 + alpha ) * real ( n * n, kind = rk ) ) )

      xval = xval + r1 * r2 * r3 * ( xval - x(i-2) )

    end if

    call gegenbauer_ss_root ( xval, n, dp2, p1, c )

    x(i) = xval
    w(i) = cc / ( dp2 * p1 )

  end do
!
!  Reverse the data.
!
  x(1:n) = x(n:1:-1)
  w(1:n) = w(n:1:-1)

  return
end
subroutine gegenbauer_ss_recur ( p2, dp2, p1, x, n, c )

!*****************************************************************************80
!
!! gegenbauer_ss_recur(): value and derivative of a Gegenbauer polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Input:
!
!    real ( kind = rk ) X, the point at which polynomials are evaluated.
!
!    integer N, the order of the polynomial.
!
!    real ( kind = rk ) C(N), the recursion coefficients.
!
!  Output:
!
!    real ( kind = rk ) P2, the value of J(N)(X).
!
!    real ( kind = rk ) DP2, the value of J'(N)(X).
!
!    real ( kind = rk ) P1, the value of J(N-1)(X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) c(n)
  real ( kind = rk ) dp0
  real ( kind = rk ) dp1
  real ( kind = rk ) dp2
  integer i
  real ( kind = rk ) p0
  real ( kind = rk ) p1
  real ( kind = rk ) p2
  real ( kind = rk ) x

  p1 = 1.0D+00
  dp1 = 0.0D+00

  p2 = x
  dp2 = 1.0D+00

  do i = 2, n

    p0 = p1
    dp0 = dp1

    p1 = p2
    dp1 = dp2

    p2 = x * p1 - c(i) * p0
    dp2 = x * dp1 + p1 - c(i) * dp0

  end do

  return
end
subroutine gegenbauer_ss_root ( x, n, dp2, p1, c )

!*****************************************************************************80
!
!! gegenbauer_ss_root() improves an approximate root of a Gegenbauer polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Input:
!
!    real ( kind = rk ) X, the approximate root.
!
!    integer N, the order of the polynomial.
!
!    real ( kind = rk ) C(N), the recursion coefficients.
!
!  Output:
!
!    real ( kind = rk ) X, the improved root.
!
!    real ( kind = rk ) DP2, the value of J'(N)(X).
!
!    real ( kind = rk ) P1, the value of J(N-1)(X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) c(n)
  real ( kind = rk ) d
  real ( kind = rk ) dp2
  real ( kind = rk ) eps
  real ( kind = rk ) p1
  real ( kind = rk ) p2
  integer step
  integer, parameter :: step_max = 10
  real ( kind = rk ) x

  eps = epsilon ( eps )

  do step = 1, step_max

    call gegenbauer_ss_recur ( p2, dp2, p1, x, n, c )

    d = p2 / dp2
    x = x - d

    if ( abs ( d ) <= eps * ( abs ( x ) + 1.0D+00 ) ) then
      return
    end if

  end do

  return
end
subroutine gegenbauer_to_monomial_matrix ( n, alpha, A )

!*****************************************************************************80
!
!! gegenbauer_to_monomial_matrix(): Gegenbauer to monomial conversion matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N: the order of A.
!
!    real alpha: the parameter.
!
!  Output:
!
!    real ( kind = rk ) A(N,N): the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) A(n,n)
  real ( kind = rk ) alpha
  real ( kind = rk ) c1
  real ( kind = rk ) c2
  integer i
  integer j
  integer nn

  if ( n <= 0 ) then
    return
  end if

  do j = 1, n
    do i = 1, n
      A(i,j) = 0.0D+00
    end do
  end do

  A(1,1) = 1.0D+00

  if ( n == 1 ) then
    return
  end if

  A(2,2) = 2.0D+00 * alpha
!
!  Perform convex sum.
!  Translating "(n+1) C(n+1) = 2 (n+alpha) x C(n) - ( n + 2 alpha - 1 ) C(n-1)"
!  drove me nuts, between indexing at 1 rather than 0, and dealing with
!  the interpretation of "n+1", because I now face the rare "off by 2" error!
!
  do j = 3, n
    nn = j - 2
    c1 = ( 2.0 * nn + 2.0 * alpha       ) / real ( nn + 1, kind = rk )
    c2 = (     - nn - 2.0 * alpha + 1.0 ) / real ( nn + 1, kind = rk )
    A(2:j,j)   =              c1 * A(1:j-1,j-1)
    A(1:j-2,j) = A(1:j-2,j) + c2 * A(1:j-2,j-2)
  end do

  return
end
subroutine hyper_2f1_values ( n_data, a, b, c, x, fx )

!*****************************************************************************80
!
!! hyper_2f1_values() returns some values of the hypergeometric function 2F1.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      fx = Hypergeometric2F1 [ a, b, c, x ]
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 September 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Input:
!
!    integer N_DATA.  The user sets N_DATA to 0 before the first call. 
!
!  Output:
!
!    integer N_DATA.  The routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    real ( kind = rk ) A, B, C, X, the parameters.
!
!    real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 24

  real ( kind = rk ) a
  real ( kind = rk ), save, dimension ( n_max ) :: a_vec = (/ &
   -2.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    2.5D+00, &
   -2.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    2.5D+00, &
   -2.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    2.5D+00, &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00, &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00, &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00 /)
  real ( kind = rk ) b
  real ( kind = rk ), save, dimension ( n_max ) :: b_vec = (/ &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00, &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00, &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00 /)
  real ( kind = rk ) c
  real ( kind = rk ), save, dimension ( n_max ) :: c_vec = (/ &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
   -5.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    4.5D+00, &
   -5.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    4.5D+00, &
   -5.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    4.5D+00 /)
  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.72356129348997784913D+00, &
    0.97911109345277961340D+00, &
    1.0216578140088564160D+00, &
    1.4051563200112126405D+00, &
    0.46961431639821611095D+00, &
    0.95296194977446325454D+00, &
    1.0512814213947987916D+00, &
    2.3999062904777858999D+00, &
    0.29106095928414718320D+00, &
    0.92536967910373175753D+00, &
    1.0865504094806997287D+00, &
    5.7381565526189046578D+00, &
    15090.669748704606754D+00, &
   -104.31170067364349677D+00, &
    21.175050707768812938D+00, &
    4.1946915819031922850D+00, &
    1.0170777974048815592D+10, &
   -24708.635322489155868D+00, &
    1372.2304548384989560D+00, &
    58.092728706394652211D+00, &
    5.8682087615124176162D+18, &
   -4.4635010147295996680D+08, &
    5.3835057561295731310D+06, &
    20396.913776019659426D+00 /)
  integer n_data
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( n_max ) :: x_vec = (/ &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.55D+00, &
    0.55D+00, &
    0.55D+00, &
    0.55D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.55D+00, &
    0.55D+00, &
    0.55D+00, &
    0.55D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    c = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    c = c_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine imtqlx ( n, d, e, z )

!*****************************************************************************80
!
!! imtqlx() diagonalizes a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This routine is a slightly modified version of the EISPACK routine to
!    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
!
!    The authors thank the authors of EISPACK for permission to use this
!    routine.
!
!    It has been modified to produce the product Q' * Z, where Z is an input
!    vector and Q is the orthogonal matrix diagonalizing the input matrix.
!    The changes consist (essentially) of applying the orthogonal
!    transformations directly to Z as they are generated.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    This version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!    Roger Martin, James Wilkinson,
!    The Implicit QL Algorithm,
!    Numerische Mathematik,
!    Volume 12, Number 5, December 1968, pages 377-383.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real ( kind = rk ) D(N), the diagonal entries of the matrix.
!    On output, the information in D has been overwritten.
!
!    Input/output, real ( kind = rk ) E(N), the subdiagonal entries of the
!    matrix, in entries E(1) through E(N-1).  On output, the information in
!    E has been overwritten.
!
!    Input/output, real ( kind = rk ) Z(N).  On input, a vector.  On output,
!    the value of Q' * Z, where Q is the matrix that diagonalizes the
!    input symmetric tridiagonal matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d(n)
  real ( kind = rk ) e(n)
  real ( kind = rk ) f
  real ( kind = rk ) g
  integer i
  integer ii
  integer, parameter :: itn = 30
  integer j
  integer k
  integer l
  integer m
  integer mml
  real ( kind = rk ) p
  real ( kind = rk ) prec
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) z(n)

  prec = epsilon ( prec )

  if ( n == 1 ) then
    return
  end if

  e(n) = 0.0D+00

  do l = 1, n

    j = 0

    do

      do m = l, n

        if ( m == n ) then
          exit
        end if

        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
          exit
        end if

      end do

      p = d(l)

      if ( m == l ) then
        exit
      end if

      if ( itn <= j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'imtqlx(): Fatal error!'
        write ( *, '(a)' ) '  Iteration limit exceeded.'
        write ( *, '(a,i8)' ) '  J = ', j
        write ( *, '(a,i8)' ) '  L = ', l
        write ( *, '(a,i8)' ) '  M = ', m
        write ( *, '(a,i8)' ) '  N = ', n
        stop 1
      end if

      j = j + 1
      g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
      r =  sqrt ( g * g + 1.0D+00 )
      g = d(m) - p + e(l) / ( g + sign ( r, g ) )
      s = 1.0D+00
      c = 1.0D+00
      p = 0.0D+00
      mml = m - l

      do ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)

        if ( abs ( g ) <= abs ( f ) ) then
          c = g / f
          r =  sqrt ( c * c + 1.0D+00 )
          e(i+1) = f * r
          s = 1.0D+00 / r
          c = c * s
        else
          s = f / g
          r =  sqrt ( s * s + 1.0D+00 )
          e(i+1) = g * r
          c = 1.0D+00 / r
          s = s * c
        end if

        g = d(i+1) - p
        r = ( d(i) - g ) * s + 2.0D+00 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
        f = z(i+1)
        z(i+1) = s * z(i) + c * f
        z(i) = c * z(i) - s * f

      end do

      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.0D+00

    end do

  end do
!
!  Sorting.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n
      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if
    end do

    if ( k /= i ) then
      d(k) = d(i)
      d(i) = p
      p = z(i)
      z(i) = z(k)
      z(k) = p
    end if

  end do

  return
end
subroutine monomial_to_gegenbauer_matrix ( n, alpha, A )

!*****************************************************************************80
!
!! monomial_to_gegenbauer_matrix(): monomial to Gegenbauer conversion matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n: the order of A.
!
!    real ( kind = rk ) alpha: the parameter.
!
!  Output:
!
!    real ( kind = rk ) A[N,N]: the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) A(0:n-1,0:n-1)
  real ( kind = rk ) alpha
  real ( kind = rk ) bot
  integer i
  integer ilo
  integer j
  real ( kind = rk ) r8_factorial
  real ( kind = rk ) top

  A(0:n-1,0:n-1) = 0.0D+00

  if ( n <= 0 ) then
    return
  end if

  do j = 0, n - 1

    ilo = mod ( j, 2 )

    do i = ilo, j, 2

      top = ( i + alpha ) * r8_factorial ( j ) * gamma ( alpha )

      bot = 2.0D+00 ** j * r8_factorial ( ( j - i ) / 2 ) &
        * gamma ( ( j + i ) / 2 + alpha + 1.0D+00 )

      A(i,j) = top / bot

    end do
  end do

  return
end
subroutine psi_values ( n_data, x, fx )

!*****************************************************************************80
!
!! psi_values() returns some values of the Psi or Digamma function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      PolyGamma[x]
!
!    or
!
!      PolyGamma[0,x]
!
!    PSI(X) = d ln ( Gamma ( X ) ) / d X = Gamma'(X) / Gamma(X)
!
!    PSI(1) = -Euler's constant.
!
!    PSI(X+1) = PSI(X) + 1 / X.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 June 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Input:
!
!    integer N_DATA.  The user sets N_DATA to 0 before the first call.  
!
!  Output:
!
!    integer N_DATA.  The routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    real ( kind = rk ) X, the argument of the function.
!
!    real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 20

  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
   -10.42375494041108D+00, &
    -5.289039896592188D+00, &
    -3.502524222200133D+00, &
    -2.561384544585116D+00, &
    -1.963510026021423D+00, &
    -1.540619213893190D+00, &
    -1.220023553697935D+00, &
    -0.9650085667061385D+00, &
    -0.7549269499470514D+00, &
    -0.5772156649015329D+00, &
    -0.4237549404110768D+00, &
    -0.2890398965921883D+00, &
    -0.1691908888667997D+00, &
    -0.6138454458511615D-01, &
     0.3648997397857652D-01, &
     0.1260474527734763D+00, &
     0.2085478748734940D+00, &
     0.2849914332938615D+00, &
     0.3561841611640597D+00, &
     0.4227843350984671D+00 /)
  integer n_data
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1D+00, &
    0.2D+00, &
    0.3D+00, &
    0.4D+00, &
    0.5D+00, &
    0.6D+00, &
    0.7D+00, &
    0.8D+00, &
    0.9D+00, &
    1.0D+00, &
    1.1D+00, &
    1.2D+00, &
    1.3D+00, &
    1.4D+00, &
    1.5D+00, &
    1.6D+00, &
    1.7D+00, &
    1.8D+00, &
    1.9D+00, &
    2.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

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
!  Input:
!
!    integer N: the argument of the factorial function.
!    If N is less than 1, the function value is returned as 1.
!
!  Output:
!
!    real ( kind = rk ) R8_FACTORIAL: the factorial of N.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_factorial
  integer i
  integer n
  real ( kind = rk ) value

  value = 1.0D+00

  do i = 1, n
    value = value * real ( i, kind = rk )
  end do

  r8_factorial = value

  return
end
subroutine r8_hyper_2f1 ( a_input, b_input, c_input, x_input, hf )

!*****************************************************************************80
!
!! r8_hyper_2f1() evaluates the hypergeometric function F(A,B,C,X).
!
!  Discussion:
!
!    A minor bug was corrected.  The HW variable, used in several places as
!    the "old" value of a quantity being iteratively improved, was not
!    being initialized.  JVB, 11 February 2008.
!
!    The original version of this program allowed the input arguments to
!    be modified, although they were restored to their input values before exit.
!    This is unacceptable if the input arguments are allowed to be constants.
!    The code has been modified so that the input arguments are never modified.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 October 2008
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    This version by John Burkardt.
!
!    The original FORTRAN77 version of this routine is copyrighted by
!    Shanjie Zhang and Jianming Jin.  However, they give permission to
!    incorporate this routine into a user program provided that the copyright
!    is acknowledged.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Input:
!
!    real ( kind = rk ) A_INPUT, B_INPUT, C_INPUT, X_INPUT,
!    the arguments of the function.  The user is allowed to pass these
!    values as constants or variables.
!    C_INPUT must not be equal to a nonpositive integer.
!    X_INPUT < 1.
!
!  Output:
!
!    real ( kind = rk ) HF, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) a_input
  real ( kind = rk ) a0
  real ( kind = rk ) aa
  real ( kind = rk ) b
  real ( kind = rk ) b_input
  real ( kind = rk ) bb
  real ( kind = rk ) c
  real ( kind = rk ) c_input
  real ( kind = rk ) c0
  real ( kind = rk ) c1
  real ( kind = rk ), parameter :: el = 0.5772156649015329D+00
  real ( kind = rk ) eps
  real ( kind = rk ) f0
  real ( kind = rk ) f1
  real ( kind = rk ) g0
  real ( kind = rk ) g1
  real ( kind = rk ) g2
  real ( kind = rk ) g3
  real ( kind = rk ) ga
  real ( kind = rk ) gabc
  real ( kind = rk ) gam
  real ( kind = rk ) gb
  real ( kind = rk ) gbm
  real ( kind = rk ) gc
  real ( kind = rk ) gca
  real ( kind = rk ) gcab
  real ( kind = rk ) gcb
  real ( kind = rk ) gm
  real ( kind = rk ) hf
  real ( kind = rk ) hw
  integer j
  integer k
  logical l0
  logical l1
  logical l2
  logical l3
  logical l4
  logical l5
  integer m
  integer nm
  real ( kind = rk ) pa
  real ( kind = rk ) pb
  real ( kind = rk ) r
  real ( kind = rk ) r0
  real ( kind = rk ) r1
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) r8_psi
  real ( kind = rk ) rm
  real ( kind = rk ) rp
  real ( kind = rk ) sm
  real ( kind = rk ) sp
  real ( kind = rk ) sp0
  real ( kind = rk ) x
  real ( kind = rk ) x_input
  real ( kind = rk ) x1
!
!  Immediately copy the input arguments!
!
  a = a_input
  b = b_input
  c = c_input
  x = x_input

  l0 = ( c == aint ( c ) ) .and. ( c < 0.0D+00 )
  l1 = ( 1.0D+00 - x < 1.0D-15 ) .and. ( c - a - b <= 0.0D+00 )
  l2 = ( a == aint ( a ) ) .and. ( a < 0.0D+00 )
  l3 = ( b == aint ( b ) ) .and. ( b < 0.0D+00 )
  l4 = ( c - a == aint ( c - a ) ) .and. ( c - a <= 0.0D+00 )
  l5 = ( c - b == aint ( c - b ) ) .and. ( c - b <= 0.0D+00 )

  if ( l0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'r8_hyper_2f1(): Fatal error!'
    write ( *, '(a)' ) '  Integral C < 0.'
    write ( *, '(a)' ) '  The hypergeometric series is divergent.'
    stop 1
  end if

  if ( l1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'r8_hyper_2f1(): Fatal error!'
    write ( *, '(a)' ) '  The hypergeometric series is divergent.'
    write ( *, '(a)' ) '  1 - X < 0, C - A - B < 0.'
    stop 1
  end if

  if ( 0.95D+00 < x ) then
    eps = 1.0D-08
  else
    eps = 1.0D-15
  end if

  if ( x == 0.0D+00 .or. a == 0.0D+00 .or. b == 0.0D+00 ) then

    hf = 1.0D+00
    return

  else if ( 1.0D+00 - x == eps .and. 0.0D+00 < c - a - b ) then

    gc = gamma ( c )
    gcab = gamma ( c - a - b )
    gca = gamma ( c - a )
    gcb = gamma ( c - b )
    hf = gc * gcab / ( gca * gcb )
    return

  else if ( 1.0D+00 + x <= eps .and. abs ( c - a + b - 1.0D+00 ) <= eps ) then

    g0 = sqrt ( r8_pi ) * 2.0D+00 ** ( - a )
    g1 = gamma ( c )
    g2 = gamma ( 1.0D+00 + a / 2.0D+00 - b )
    g3 = gamma ( 0.5D+00 + 0.5D+00 * a )
    hf = g0 * g1 / ( g2 * g3 )
    return

  else if ( l2 .or. l3 ) then

    if ( l2 ) then
      nm = int ( abs ( a ) )
    end if

    if ( l3 ) then
      nm = int ( abs ( b ) )
    end if

    hf = 1.0D+00
    r = 1.0D+00

    do k = 1, nm
      r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
        / ( k * ( c + k - 1.0D+00 ) ) * x
      hf = hf + r
    end do

    return

  else if ( l4 .or. l5 ) then

    if ( l4 ) then
      nm = int ( abs ( c - a ) )
    end if

    if ( l5 ) then
      nm = int ( abs ( c - b ) )
    end if

    hf = 1.0D+00
    r  = 1.0D+00
    do k = 1, nm
      r = r * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
        / ( k * ( c + k - 1.0D+00 ) ) * x
      hf = hf + r
    end do
    hf = ( 1.0D+00 - x )**( c - a - b ) * hf
    return

  end if

  aa = a
  bb = b
  x1 = x

  if ( x < 0.0D+00 ) then
    x = x / ( x - 1.0D+00 )
    if ( a < c .and. b < a .and. 0.0D+00 < b ) then
      a = bb
      b = aa
    end if
    b = c - b
  end if

  if ( 0.75D+00 <= x ) then

    gm = 0.0D+00

    if ( abs ( c - a - b - aint ( c - a - b ) ) < 1.0D-15 ) then

      m = int ( c - a - b )
      ga = gamma ( a )
      gb = gamma ( b )
      gc = gamma ( c )
      gam = gamma ( a + m )
      gbm = gamma ( b + m )

      pa = r8_psi ( a )
      pb = r8_psi ( b )

      if ( m /= 0 ) then
        gm = 1.0D+00
      end if

      do j = 1, abs ( m ) - 1
        gm = gm * j
      end do

      rm = 1.0D+00
      do j = 1, abs ( m )
        rm = rm * j
      end do

      f0 = 1.0D+00
      r0 = 1.0D+00
      r1 = 1.0D+00
      sp0 = 0.0D+00
      sp = 0.0D+00

      if ( 0 <= m ) then

        c0 = gm * gc / ( gam * gbm )
        c1 = - gc * ( x - 1.0D+00 )**m / ( ga * gb * rm )

        do k = 1, m - 1
          r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
            / ( k * ( k - m ) ) * ( 1.0D+00 - x )
          f0 = f0 + r0
        end do

        do k = 1, m
          sp0 = sp0 + 1.0D+00 / ( a + k - 1.0D+00 ) &
            + 1.0D+00 / ( b + k - 1.0D+00 ) - 1.0D+00 / real ( k, kind = rk )
        end do

        f1 = pa + pb + sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
        hw = f1

        do k = 1, 250

          sp = sp + ( 1.0D+00 - a ) / ( k * ( a + k - 1.0D+00 ) ) &
            + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

          sm = 0.0D+00
          do j = 1, m
            sm = sm + ( 1.0D+00 - a ) &
              / ( ( j + k ) * ( a + j + k - 1.0D+00 ) ) &
              + 1.0D+00 / ( b + j + k - 1.0D+00 )
          end do

          rp = pa + pb + 2.0D+00 * el + sp + sm + log ( 1.0D+00 - x )

          r1 = r1 * ( a + m + k - 1.0D+00 ) * ( b + m + k - 1.0D+00 ) &
            / ( k * ( m + k ) ) * ( 1.0D+00 - x )

          f1 = f1 + r1 * rp

          if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
            exit
          end if

          hw = f1

        end do

        hf = f0 * c0 + f1 * c1

      else if ( m < 0 ) then

        m = - m
        c0 = gm * gc / ( ga * gb * ( 1.0D+00 - x )**m )
        c1 = - ( - 1 )**m * gc / ( gam * gbm * rm )

        do k = 1, m - 1
          r0 = r0 * ( a - m + k - 1.0D+00 ) * ( b - m + k - 1.0D+00 ) &
            / ( k * ( k - m ) ) * ( 1.0D+00 - x )
          f0 = f0 + r0
        end do

        do k = 1, m
          sp0 = sp0 + 1.0D+00 / real ( k, kind = rk )
        end do

        f1 = pa + pb - sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
        hw = f1

        do k = 1, 250

          sp = sp + ( 1.0D+00 - a ) &
            / ( k * ( a + k - 1.0D+00 ) ) &
            + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

          sm = 0.0D+00
          do j = 1, m
            sm = sm + 1.0D+00 / real ( j + k, kind = rk )
          end do

          rp = pa + pb + 2.0D+00 * el + sp - sm + log ( 1.0D+00 - x )

          r1 = r1 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
            / ( k * ( m + k ) ) * ( 1.0D+00 - x )

          f1 = f1 + r1 * rp

          if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
            exit
          end if

          hw = f1

        end do

        hf = f0 * c0 + f1 * c1

      end if

    else

      ga = gamma ( a )
      gb = gamma ( b )
      gc = gamma ( c )
      gca = gamma ( c - a )
      gcb = gamma ( c - b )
      gcab = gamma ( c - a - b )
      gabc = gamma ( a + b - c )
      c0 = gc * gcab / ( gca * gcb )
      c1 = gc * gabc / ( ga * gb ) * ( 1.0D+00 - x )**( c - a - b )
      hf = 0.0D+00
      hw = hf
      r0 = c0
      r1 = c1

      do k = 1, 250

        r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
          / ( k * ( a + b - c + k ) ) * ( 1.0D+00 - x )

        r1 = r1 * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
          / ( k * ( c - a - b + k ) ) * ( 1.0D+00 - x )

        hf = hf + r0 + r1

        if ( abs ( hf - hw ) < abs ( hf ) * eps ) then
          exit
        end if

        hw = hf

      end do

      hf = hf + c0 + c1

    end if

  else

    a0 = 1.0D+00

    if ( a < c .and. c < 2.0D+00 * a .and. b < c .and. c < 2.0D+00 * b ) then

      a0 = ( 1.0D+00 - x )**( c - a - b )
      a = c - a
      b = c - b

    end if

    hf = 1.0D+00
    hw = hf
    r = 1.0D+00

    do k = 1, 250

      r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
        / ( k * ( c + k - 1.0D+00 ) ) * x

      hf = hf + r

      if ( abs ( hf - hw ) <= abs ( hf ) * eps ) then
        exit
      end if

      hw = hf

    end do

    hf = a0 * hf

  end if

  if ( x1 < 0.0D+00 ) then
    x = x1
    c0 = 1.0D+00 / ( 1.0D+00 - x ) ** aa
    hf = c0 * hf
  end if

  a = aa
  b = bb

  if ( 120 < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'r8_hyper_2f1(): Warning!'
    write ( *, '(a)' ) '  A large number of iterations were needed.'
    write ( *, '(a)' ) '  The accuracy of the results should be checked.'
  end if

  return
end
function r8_psi ( xx )

!*****************************************************************************80
!
!! r8_psi() evaluates the function Psi(X).
!
!  Discussion:
!
!    This routine evaluates the logarithmic derivative of the
!    Gamma function,
!
!      PSI(X) = d/dx ( GAMMA(X) ) / GAMMA(X)
!             = d/dx LN ( GAMMA(X) )
!
!    for real X, where either
!
!      - XMAX1 < X < - XMIN, and X is not a negative integer,
!
!    or
!
!      XMIN < X.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    This version by John Burkardt.
!
!  Reference:
!
!    William Cody, Anthony Strecok, Henry Thacher,
!    Chebyshev Approximations for the Psi Function,
!    Mathematics of Computation,
!    Volume 27, Number 121, January 1973, pages 123-127.
!
!  Input:
!
!    real ( kind = rk ) XX, the argument of the function.
!
!  Output:
!
!    real ( kind = rk ) R8_PSI, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) aug
  real ( kind = rk ) den
  real ( kind = rk ), parameter :: four = 4.0D+00
  real ( kind = rk ), parameter :: fourth = 0.25D+00
  real ( kind = rk ), parameter :: half = 0.5D+00
  integer i
  integer n
  integer nq
  real ( kind = rk ), parameter :: one = 1.0D+00
  real ( kind = rk ), dimension ( 9 ) :: p1 = (/ &
   4.5104681245762934160D-03, &
   5.4932855833000385356D+00, &
   3.7646693175929276856D+02, &
   7.9525490849151998065D+03, &
   7.1451595818951933210D+04, &
   3.0655976301987365674D+05, &
   6.3606997788964458797D+05, &
   5.8041312783537569993D+05, &
   1.6585695029761022321D+05 /)
  real ( kind = rk ), dimension ( 7 ) :: p2 = (/ &
  -2.7103228277757834192D+00, &
  -1.5166271776896121383D+01, &
  -1.9784554148719218667D+01, &
  -8.8100958828312219821D+00, &
  -1.4479614616899842986D+00, &
  -7.3689600332394549911D-02, &
  -6.5135387732718171306D-21 /)
  real ( kind = rk ), parameter :: piov4 = 0.78539816339744830962D+00
  real ( kind = rk ), dimension ( 8 ) :: q1 = (/ &
   9.6141654774222358525D+01, &
   2.6287715790581193330D+03, &
   2.9862497022250277920D+04, &
   1.6206566091533671639D+05, &
   4.3487880712768329037D+05, &
   5.4256384537269993733D+05, &
   2.4242185002017985252D+05, &
   6.4155223783576225996D-08 /)
  real ( kind = rk ), dimension ( 6 ) :: q2 = (/ &
   4.4992760373789365846D+01, &
   2.0240955312679931159D+02, &
   2.4736979003315290057D+02, &
   1.0742543875702278326D+02, &
   1.7463965060678569906D+01, &
   8.8427520398873480342D-01 /)
  real ( kind = rk ) r8_psi
  real ( kind = rk ) sgn
  real ( kind = rk ), parameter :: three = 3.0D+00
  real ( kind = rk ) upper
  real ( kind = rk ) w
  real ( kind = rk ) x
  real ( kind = rk ), parameter :: x01 = 187.0D+00
  real ( kind = rk ), parameter :: x01d = 128.0D+00
  real ( kind = rk ), parameter :: x02 = 6.9464496836234126266D-04
  real ( kind = rk ), parameter :: xinf = 1.70D+38
  real ( kind = rk ), parameter :: xlarge = 2.04D+15
  real ( kind = rk ), parameter :: xmax1 = 3.60D+16
  real ( kind = rk ), parameter :: xmin1 = 5.89D-39
  real ( kind = rk ), parameter :: xsmall = 2.05D-09
  real ( kind = rk ) xx
  real ( kind = rk ) z
  real ( kind = rk ), parameter :: zero = 0.0D+00

  x = xx
  w = abs ( x )
  aug = zero
!
!  Check for valid arguments, then branch to appropriate algorithm.
!
  if ( xmax1 <= - x .or. w < xmin1 ) then

    if ( zero < x ) then
      r8_psi = - xinf
    else
      r8_psi = xinf
    end if

    return
  end if

  if ( x < half ) then
!
!  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
!  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
!
    if ( w <= xsmall ) then

      aug = - one / x
!
!  Argument reduction for cotangent.
!
    else

      if ( x < zero ) then
        sgn = piov4
      else
        sgn = - piov4
      end if

      w = w - real ( int ( w ), kind = rk )
      nq = int ( w * four )
      w = four * ( w - real ( nq, kind = rk ) * fourth )
!
!  W is now related to the fractional part of 4.0 * X.
!  Adjust argument to correspond to values in the first
!  quadrant and determine the sign.
!
      n = nq / 2

      if ( n + n /= nq ) then
        w = one - w
      end if

      z = piov4 * w

      if ( mod ( n, 2 ) /= 0 ) then
        sgn = - sgn
      end if
!
!  Determine the final value for  -pi * cotan(pi*x).
!
      n = ( nq + 1 ) / 2
      if ( mod ( n, 2 ) == 0 ) then
!
!  Check for singularity.
!
        if ( z == zero ) then

          if ( zero < x ) then
            r8_psi = -xinf
          else
            r8_psi = xinf
          end if

          return
        end if

        aug = sgn * ( four / tan ( z ) )

      else

        aug = sgn * ( four * tan ( z ) )

      end if

    end if

    x = one - x

  end if
!
!  0.5 <= X <= 3.0.
!
  if ( x <= three ) then

    den = x
    upper = p1(1) * x
    do i = 1, 7
      den = ( den + q1(i) ) * x
      upper = ( upper + p1(i+1) ) * x
    end do
    den = ( upper + p1(9) ) / ( den + q1(8) )
    x = ( x - x01 / x01d ) - x02
    r8_psi = den * x + aug
    return

  end if
!
!  3.0 < X.
!
  if ( x < xlarge ) then
    w = one / ( x * x )
    den = w
    upper = p2(1) * w
    do i = 1, 5
      den = ( den + q2(i) ) * w
      upper = ( upper + p2(i+1) ) * w
    end do
    aug = ( upper + p2(7) ) / ( den + q2(6) ) - half / x + aug
  end if

  r8_psi = aug + log ( x )

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! r8vec_print() prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real ( kind = rk ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end

