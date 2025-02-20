function digamma ( x, ifault )

!*****************************************************************************80
!
!! digamma() calculates DIGAMMA ( X ) = d ( LOG ( GAMMA ( X ) ) ) / dX
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 August 2021
!
!  Author:
!
!    Original FORTRAN77 version by Jose Bernardo.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jose Bernardo,
!    Algorithm AS 103:
!    Psi ( Digamma ) Function,
!    Applied Statistics,
!    Volume 25, Number 3, 1976, pages 315-317.
!
!  Input:
!
!    real ( kind = rk ) X, the argument of the digamma function.
!    0 < X.
!
!  Output:
!
!    integer IFAULT, error flag.
!    0, no error.
!    1, X <= 0.
!
!    real ( kind = rk ) DIGAMMA, the value of the digamma function at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: c = 8.5D+00
  real ( kind = rk ), parameter :: euler_mascheroni = 0.57721566490153286060D+00
  real ( kind = rk ) digamma
  integer ifault
  real ( kind = rk ) r
  real ( kind = rk ) x
  real ( kind = rk ) x2
!
!  Check the input.
!
  if ( x <= 0.0D+00 ) then
    digamma = 0.0D+00
    ifault = 1
    return
  end if
!
!  Initialize.
!
  ifault = 0
!
!  Approximation for small argument.
!
  if ( x <= 0.000001D+00 ) then
    digamma = - euler_mascheroni - 1.0D+00 / x + 1.6449340668482264365D+00 * x
    return
  end if
!
!  Reduce to DIGAMA(X + N).
!
  digamma = 0.0D+00
  x2 = x

  do while ( x2 < c )
    digamma = digamma - 1.0D+00 / x2
    x2 = x2 + 1.0D+00
  end do
!
!  Use Stirling's (actually de Moivre's) expansion.
!
  r = 1.0D+00 / x2

  digamma = digamma + log ( x2 ) - 0.5D+00 * r

  r = r * r

  digamma = digamma &
    - r * ( 1.0D+00 / 12.0D+00 &
    - r * ( 1.0D+00 / 120.0D+00 &
    - r * ( 1.0D+00 / 252.0D+00 &
    - r * ( 1.0D+00 / 240.0D+00 &
    - r * ( 1.0D+00 / 132.0D+00 ) ) ) ) )

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
!    28 August 2021
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
!    integer N_DATA.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    real ( kind = rk ) X, the argument of the function.
!
!    real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 11

  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
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

