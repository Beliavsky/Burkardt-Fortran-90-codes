function local_min ( a, b, eps, t, f, x, calls )

!*****************************************************************************80
!
!! local_min() seeks a local minimum of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    If the function F is defined on the interval (A,B), then local_min
!    finds an approximation X to the point at which F attatains its minimum
!    (or the appropriate limit point), and returns the value of F at X.
!
!    T and EPS define a tolerance TOL = EPS * abs ( X ) + T.
!    F is never evaluated at two points closer than TOL.  
!
!    If F is delta-unimodal for some delta less than TOL, the X approximates
!    the global minimum of F with an error less than 3*TOL.
!
!    If F is not delta-unimodal, then X may approximate a local, but 
!    perhaps non-global, minimum.
!
!    The method used is a combination of golden section search and
!    successive parabolic interpolation.  Convergence is never much slower
!    than that for a Fibonacci search.  If F has a continuous second
!    derivative which is positive at the minimum (which is not at A or
!    B), then, ignoring rounding errors, convergence is superlinear, and 
!    usually of the order of about 1.3247.
!
!    Thanks to Jonathan Eggleston for pointing out a correction to the 
!    golden section step, 01 July 2013.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 May 2021
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Input:
!
!    real ( kind = rk ) A, B, the endpoints of the interval.
!
!    real ( kind = rk ) EPS, a positive relative error tolerance.
!    EPS should be no smaller than twice the relative machine precision,
!    and preferably not much less than the square root of the relative
!    machine precision.
!
!    real ( kind = rk ) T, a positive absolute error tolerance.
!
!    external real ( kind = rk ) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose local minimum is being sought.
!
!  Output:
!
!    real ( kind = rk ) X, the estimated value of an abscissa
!    for which F attains a local minimum value in [A,B].
!
!    integer ( kind = 4 ) CALLS: the number of calls to F.
!
!    real ( kind = rk ) LOCAL_MIN, the value F(X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  integer ( kind = 4 ) calls
  real ( kind = rk ) d
  real ( kind = rk ) e
  real ( kind = rk ) eps
  real ( kind = rk ) f
  real ( kind = rk ) fu
  real ( kind = rk ) fv
  real ( kind = rk ) fw
  real ( kind = rk ) fx
  real ( kind = rk ) local_min
  real ( kind = rk ) m
  real ( kind = rk ) p
  real ( kind = rk ) q
  real ( kind = rk ) r
  real ( kind = rk ) sa
  real ( kind = rk ) sb
  real ( kind = rk ) t
  real ( kind = rk ) t2
  real ( kind = rk ) tol
  real ( kind = rk ) u
  real ( kind = rk ) v
  real ( kind = rk ) w
  real ( kind = rk ) x

  calls = 0
!
!  C is the square of the inverse of the golden ratio.
!
  c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )

  sa = a
  sb = b
  x = sa + c * ( b - a )
  w = x
  v = w
  e = 0.0D+00
  fx = f ( x )
  calls = calls + 1
  fw = fx
  fv = fw

  do

    m = 0.5D+00 * ( sa + sb ) 
    tol = eps * abs ( x ) + t
    t2 = 2.0D+00 * tol
!
!  Check the stopping criterion.
!
    if ( abs ( x - m ) <= t2 - 0.5D+00 * ( sb - sa ) ) then
      exit
    end if
!
!  Fit a parabola.
!
    r = 0.0D+00
    q = r
    p = q

    if ( tol < abs ( e ) ) then

      r = ( x - w ) * ( fx - fv )
      q = ( x - v ) * ( fx - fw )
      p = ( x - v ) * q - ( x - w ) * r
      q = 2.0D+00 * ( q - r )

      if ( 0.0D+00 < q ) then
        p = - p
      end if

      q = abs ( q )

      r = e
      e = d

    end if

    if ( abs ( p ) < abs ( 0.5D+00 * q * r ) .and. &
         q * ( sa - x ) < p .and. &
         p < q * ( sb - x ) ) then
!
!  Take the parabolic interpolation step.
!
      d = p / q
      u = x + d
!
!  F must not be evaluated too close to A or B.
!
      if ( ( u - sa ) < t2 .or. ( sb - u ) < t2 ) then

        if ( x < m ) then
          d = tol
        else
          d = - tol
        end if

      end if
!
!  A golden-section step.
!
    else

      if ( x < m ) then
        e = sb - x
      else
        e = sa - x
      end if

      d = c * e

    end if
!
!  F must not be evaluated too close to X.
!
    if ( tol <= abs ( d ) ) then
      u = x + d
    else if ( 0.0D+00 < d ) then
      u = x + tol
    else
      u = x - tol
    end if

    fu = f ( u )
    calls = calls + 1
!
!  Update A, B, V, W, and X.
!
    if ( fu <= fx ) then

      if ( u < x ) then
        sb = x
      else
        sa = x
      end if

      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu

    else

      if ( u < x ) then
        sa = u
      else
        sb = u
      end if

      if ( fu <= fw .or. w == x ) then
        v = w
        fv = fw
        w = u
        fw = fu
      else if ( fu <= fv .or. v == x .or. v == w ) then
        v = u
        fv = fu
      end if

    end if

  end do

  local_min = fx

  return
end

