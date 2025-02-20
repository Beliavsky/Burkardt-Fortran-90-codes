function zero_brent ( a, b, t, f, calls )

!*****************************************************************************80
!
!! zero_brent() seeks a root of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    The interval [A,B] must be a change of sign interval for F.
!    That is, F(A) and F(B) must be of opposite signs.  Then
!    assuming that F is continuous implies the existence of at least
!    one value C between A and B for which F(C) = 0.
!
!    The location of the zero is determined to within an accuracy
!    of 4 * EPSILON * abs ( C ) + 2 * T.
!
!    Thanks to Thomas Secretin for pointing out a transcription error in the
!    setting of the value of P, 11 February 2013.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 July 2021
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    This version by John Burkardt.
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
!    real ( kind = rk ) A, B: the endpoints of the change of 
!    sign interval.
!
!    real ( kind = rk ) T: a positive error tolerance.
!
!    external real ( kind = rk ) F: the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose zero is being sought.
!
!  Output:
!
!    integer CALLS: the number of calls to F.
!
!    real ( kind = rk ) zero_brent the estimated value of a zero of
!    the function F.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  integer calls
  real ( kind = rk ) d
  real ( kind = rk ) e
  real ( kind = rk ) f
  real ( kind = rk ) fa
  real ( kind = rk ) fb
  real ( kind = rk ) fc
  real ( kind = rk ) m
  real ( kind = rk ) p
  real ( kind = rk ) q
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) sa
  real ( kind = rk ) sb
  real ( kind = rk ) t
  real ( kind = rk ) tol
  real ( kind = rk ) zero_brent

  calls = 0
!
!  Make local copies of A and B.
!
  sa = a
  fa = f ( sa )
  calls = calls + 1

  sb = b
  fb = f ( sb )
  calls = calls + 1

  c = sa
  fc = fa
  e = sb - sa
  d = e

  do

    if ( abs ( fc ) < abs ( fb ) ) then

      sa = sb
      sb = c
      c = sa
      fa = fb
      fb = fc
      fc = fa

    end if

    tol = 2.0D+00 * epsilon ( sb ) * abs ( sb ) + t
    m = 0.5D+00 * ( c - sb )

    if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
      exit
    end if

    if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then

      e = m
      d = e

    else

      s = fb / fa

      if ( sa == c ) then

        p = 2.0D+00 * m * s
        q = 1.0D+00 - s

      else

        q = fa / fc
        r = fb / fc
        p = s * ( 2.0D+00 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
        q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

      end if

      if ( 0.0D+00 < p ) then
        q = - q
      else
        p = - p
      end if

      s = e
      e = d

      if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. &
        p < abs ( 0.5D+00 * s * q ) ) then
        d = p / q
      else
        e = m
        d = e
      end if

    end if

    sa = sb
    fa = fb

    if ( tol < abs ( d ) ) then
      sb = sb + d
    else if ( 0.0D+00 < m ) then
      sb = sb + tol
    else
      sb = sb - tol
    end if

    fb = f ( sb )
    calls = calls + 1

    if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
         ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
      c = sa
      fc = fa
      e = sb - sa
      d = e
    end if

  end do

  zero_brent = sb

  return
end

