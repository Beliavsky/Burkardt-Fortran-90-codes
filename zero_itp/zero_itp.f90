subroutine zero_itp ( f, a, b, epsi, k1, k2, n0, verbose, z, fz, calls )

!*****************************************************************************80
!
!! zero_itp() seeks a zero of a function using the ITP algorithm.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 March 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    function f ( x ): the name of the user-supplied function.
!
!    real a, b: the endpoints of the change of sign interval.
!
!    real epsi: error tolerance between exact and computed roots.
!
!    real k1: a parameter, with suggested value 0.2 / ( b - a ).
!
!    real k2: a parameter, typically set to 2.
!
!    integer n0: a parameter that can be set to 0 for difficult problems,
!    but is usually set to 1, to take more advantage of the secant method.
!
!    logical verbose: if true, requests additional printed output.
!
!  Output:
!
!    real z, fz: the estimated root and its function value.
!
!    integer ncalls: the number of function calls.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  integer calls
  real ( kind = rk ) delta
  real ( kind = rk ) epsi
  real ( kind = rk ), external :: f
  real ( kind = rk ) fz
  real ( kind = rk ) k1
  real ( kind = rk ) k2
  integer n0
  integer nh
  integer nmax
  real ( kind = rk ) r
  real ( kind = rk ) r8_log_2
  real ( kind = rk ) s
  real ( kind = rk ) sigma
  logical verbose
  real ( kind = rk ) xf
  real ( kind = rk ) xh
  real ( kind = rk ) xitp
  real ( kind = rk ) xt
  real ( kind = rk ) ya
  real ( kind = rk ) yb
  real ( kind = rk ) yitp
  reaL ( kind = rk ) z

  if ( b < a ) then
    c = a
    a = b
    b = c
  end if

  ya = f ( a )
  yb = f ( b )
  if ( 0.0D+00 < ya * yb ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'zero_itp(): Fatal error!'
    write ( *, '(a)' ) '  f(a) and f(b) have sign sign.'
    stop ( 1 )
  end if
!
!  Modify f(x) so that y(a) < 0, 0 < y(b)
!
  if ( 0.0D+00 < ya ) then
    s = -1.0D+00
    ya = - ya
    yb = - yb
  else
    s = +1.0D+00
  end if
  
  nh = ceiling ( r8_log_2 ( ( b - a ) / 2.0D+00 / epsi ) )
  nmax = nh + n0

  calls = 0

  if ( verbose ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) '  User has requested additional verbose output.'
    write ( *, '(a)' ) '  step   [a,    b]    x    f(x)'
    write ( *, '(a)' ) ''
  end if

  do while ( 2.0D+00 * epsi < ( b - a ) )
!
!  Calculate parameters.
!
    xh = 0.5D+00 * ( a + b )
    r = epsi * 2.0D+00 ** ( nmax - calls ) - 0.5D+00 * ( b - a )
    delta = k1 * ( b - a ) ** k2
!
!  Interpolation.
!
    xf = ( yb * a - ya * b ) / ( yb - ya )
!
!  Truncation.
!
    sigma = sign ( 1.0D+00, xh - xf )
    if ( delta < abs ( xh - xf ) ) then
      xt = xf + sigma * delta
    else
      xt = xh
    end if
!
!  Projection.
!
    if ( abs ( xt - xh ) <= r ) then
      xitp = xt
    else
      xitp = xh - sigma * r
    end if
!
!  Update the interval.
!
    yitp = s * f ( xitp )

    if ( verbose ) then
      write ( *, '(2x,i4,2x,f10.6,2x,f10.6,2x,f10.6,2x,g14.6)' ) &
        calls, a, b, xitp, yitp
    end if

    if ( 0.0D+00 < yitp ) then
      b = xitp
      yb = yitp
    else if ( yitp < 0.0D+00 ) then
      a = xitp
      ya = yitp
    else
      a = xitp
      b = xitp
      exit
    end if

    calls = calls + 1

  end do

  z = 0.5D+00 * ( a + b )
  fz = f ( z )

  return
end
function r8_log_2 ( x )

!*****************************************************************************80
!
!! r8_log_2() returns the logarithm base 2 of an R8.
!
!  Discussion:
!
!    value = Log ( |X| ) / Log ( 2.0 )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) X, the number whose base 2 logarithm is desired.
!    X should not be 0.
!
!  Output:
!
!    real ( kind = rk ) R8_LOG_2, the logarithm base 2 of the absolute
!    value of X.  It should be true that |X| = 2^R8_LOG_2.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_log_2
  real ( kind = rk ) x

  if ( x == 0.0D+00 ) then
    r8_log_2 = - huge ( x )
  else
    r8_log_2 = log ( abs ( x ) ) / log ( 2.0D+00 )
  end if

  return
end

