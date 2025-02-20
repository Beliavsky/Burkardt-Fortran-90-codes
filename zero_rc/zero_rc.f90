subroutine zero_rc ( a, b, t, arg, status, value )

!*****************************************************************************80
!
!! zero_rc() seeks a root of a function F(X) using reverse communication.
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
!    The routine is a revised version of the Brent zero finder 
!    algorithm, using reverse communication.
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
!    John Burkardt
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
!    real ( kind = rk ) A, B, the endpoints of the change of sign interval.
!
!    real ( kind = rk ) T, a positive error tolerance.
!
!    integer STATUS, used to communicate between the user 
!    and the routine.  The user only sets STATUS to zero on the first call, 
!    to indicate that this is a startup call.
!
!    real ( kind = rk ) VALUE, the function value at ARG, as requested
!    by the routine on the previous call.
!
!  Output:
!
!    real ( kind = rk ) ARG, the currently considered point.  For the next call,
!    the user is requested to evaluate the function at ARG, and return
!    the value in VALUE.  On return with STATUS zero, ARG is the routine's
!    estimate for the function's zero.
!
!    integer STATUS, used to communicate between the user 
!    and the routine.  The routine returns STATUS positive to request 
!    that the function be evaluated at ARG, or returns STATUS as 0, to 
!    indicate that the iteration is complete and that ARG is the estimated zero.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) arg
  real ( kind = rk ) b
  real ( kind = rk ), save :: c
  real ( kind = rk ), save :: d
  real ( kind = rk ), save :: e
  real ( kind = rk ), save :: fa
  real ( kind = rk ), save :: fb
  real ( kind = rk ), save :: fc
  real ( kind = rk ) m
  real ( kind = rk ) p
  real ( kind = rk ) q
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ), save :: sa
  real ( kind = rk ), save :: sb
  integer status
  real ( kind = rk ) t
  real ( kind = rk ) tol
  real ( kind = rk ) value
!
!  Input STATUS = 0.
!  Initialize, request F(A).
!
  if ( status == 0 ) then

    sa = a
    sb = b
    e = sb - sa
    d = e

    status = 1
    arg = a
    return
!
!  Input STATUS = 1.
!  Receive F(A), request F(B).
!
  else if ( status == 1 ) then

    fa = value

    status = 2
    arg = sb
    return
!
!  Input STATUS = 2
!  Receive F(B).
!
  else if ( status == 2 ) then

    fb = value

    if ( 0.0D+00 < fa * fb ) then
      status = -1
      return
    end if

    c = sa
    fc = fa

  else

    fb = value

    if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
         ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
      c = sa
      fc = fa
      e = sb - sa
      d = e
    end if

  end if
!
!  Compute the next point at which a function value is requested.
!
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
    status = 0
    arg = sb
    return
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

  arg = sb
  status = status + 1

  return
end

