function glomin ( a, b, c, m, machep, e, t, f, x, calls )

!*****************************************************************************80
!
!! glomin() seeks a global minimum of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    This function assumes: 
!    * F(X) is twice continuously differentiable over [A,B];
!    * F''(X) <= M for all X in [A,B];
!    * the user can supply the value of this upper bound M.
!
!    Thanks to Hans Bieshaar for supplying several corrections to the code,
!    28 May 2021
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 May 2021
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
!    It must be the case that A < B.
!
!    real ( kind = rk ) C, an initial guess for the global
!    minimizer.  If no good guess is known, C = A or B is acceptable.
!
!    real ( kind = rk ) M, the bound on the second derivative.
!
!    real ( kind = rk ) MACHEP, an estimate for the relative machine
!    precision.
!
!    real ( kind = rk ) E, a positive tolerance, a bound for the
!    absolute error in the evaluation of F(X) for any X in [A,B].
!
!    real ( kind = rk ) T, a positive error tolerance.
!
!    external real ( kind = rk ) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose global minimum is being sought.
!
!  Output:
!
!    real ( kind = rk ) X, the estimated value of the abscissa
!    for which F attains its global minimum value in [A,B].
!
!    integer CALLS, the number of function calls.
!
!    real ( kind = rk ) GLOMIN, the value F(X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) a0
  real ( kind = rk ) a2
  real ( kind = rk ) a3
  real ( kind = rk ) b
  real ( kind = rk ) c
  integer calls
  real ( kind = rk ) d0
  real ( kind = rk ) d1
  real ( kind = rk ) d2
  real ( kind = rk ) e
  real ( kind = rk ) f
  logical ( kind = 4 ) force_first
  real ( kind = rk ) glomin
  real ( kind = rk ) h
  integer k
  real ( kind = rk ) m
  real ( kind = rk ) m2
  real ( kind = rk ) machep
  real ( kind = rk ) p
  real ( kind = rk ) q
  real ( kind = rk ) qs
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) x
  real ( kind = rk ) y
  real ( kind = rk ) y0
  real ( kind = rk ) y1
  real ( kind = rk ) y2
  real ( kind = rk ) y3
  real ( kind = rk ) yb
  real ( kind = rk ) z0
  real ( kind = rk ) z1
  real ( kind = rk ) z2

  calls = 0
  a0 = b
  x = a0
  a2 = a
  y0 = f ( b )
  calls = calls + 1
  yb = y0
  y2 = f ( a )
  calls = calls + 1
  y = y2

  if ( y0 < y ) then
    y = y0
  else
    x = a
  end if

  if ( m <= 0.0D+00 .or. b <= a ) then
    glomin = y
    return
  end if

  m2 = 0.5D+00 * ( 1.0D+00 + 16.0D+00 * machep ) * m

  if ( c <= a .or. b <= c ) then
    c = 0.5D+00 * ( a + b )
  end if

  y1 = f ( c )
  calls = calls + 1
  k = 3
  d0 = a2 - c
  h = 9.0D+00 / 11.0D+00

  if ( y1 < y ) then
    x = c
    y = y1
  end if

  do

    d1 = a2 - a0
    d2 = c - a0
    z2 = b - a2
    z0 = y2 - y1
    z1 = y2 - y0
    r = d1 * d1 * z0 - d0 * d0 * z1
    p = r
    qs = 2.0D+00 * ( d0 * z1 - d1 * z0 )
    q = qs
!
!  Loop control corrected by Hans Bieshaar, 28 May 2021.
!
    force_first = .true.

    if ( 100000 < k .and. y < y2 ) then
      k = mod ( 1611 * k, 1048576 )
      q = 1.0D+00
      r = ( b - a ) * 0.00001D+00 * real ( k, kind = rk )
      force_first = .false.
    end if

    do while ( r < z2 .or. force_first )

      force_first = .false.

      if ( q * ( r * ( yb - y2 ) + z2 * q * ( ( y2 - y ) + t ) ) < &
        z2 * m2 * r * ( z2 * q - r ) ) then

        a3 = a2 + r / q
        y3 = f ( a3 )
        calls = calls + 1

        if ( y3 < y ) then
          x = a3
          y = y3
        end if

      end if

      k = mod ( 1611 * k, 1048576 )
      q = 1.0D+00
      r = ( b - a ) * 0.00001D+00 * real ( k, kind = rk )

    end do

    r = m2 * d0 * d1 * d2
    s = sqrt ( ( ( y2 - y ) + t ) / m2 )
    h = 0.5D+00 * ( 1.0D+00 + h )
    p = h * ( p + 2.0D+00 * r * s )
!
!  Correction supplied by Hans Bieshaar, 27 May 2021.
!
    q = r + 0.5D+00 * qs
    r = - 0.5D+00 * ( d0 + ( z0 + 2.01D+00 * e ) / ( d0 * m2 ) )

    if ( r < s .or. d0 < 0.0D+00 ) then
      r = a2 + s
    else
      r = a2 + r
    end if

    if ( 0.0D+00 < p * q ) then
      a3 = a2 + p / q
    else
      a3 = r
    end if

    do

      a3 = max ( a3, r )

      if ( b <= a3 ) then
        a3 = b
        y3 = yb
      else
        y3 = f ( a3 )
        calls = calls + 1
      end if

      if ( y3 < y ) then
        x = a3
        y = y3
      end if

      d0 = a3 - a2

      if ( a3 <= r ) then
        exit
      end if

      p = 2.0D+00 * ( y2 - y3 ) / ( m * d0 )

      if ( ( 1.0D+00 + 9.0D+00 * machep ) * d0 <= abs ( p ) ) then
        exit
      end if

      if ( 0.5D+00 * m2 * ( d0 * d0 + p * p ) <= &
        ( y2 - y ) + ( y3 - y ) + 2.0D+00 * t ) then
        exit
      end if

      a3 = 0.5D+00 * ( a2 + a3 )
      h = 0.9D+00 * h

    end do

    if ( b <= a3 ) then
      exit
    end if

    a0 = c
    c = a2
    a2 = a3
    y0 = y1
    y1 = y2
    y2 = y3

  end do

  glomin = y

  return
end

