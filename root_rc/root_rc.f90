function root_rc ( x, fx, ferr, xerr, q )

!*****************************************************************************80
!
!! root_rc() solves a single nonlinear equation using reverse communication.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 January 2013
!
!  Author:
!
!    Original FORTRAN77 version by Gaston Gonnet.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!     
!    Gaston Gonnet,
!    On the Structure of Zero Finders,
!    BIT Numerical Mathematics,
!    Volume 17, Number 2, June 1977, pages 170-183.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, an estimate for the root.  On the first
!    call, this must be a value chosen by the user.  Thereafter, it may
!    be a value chosen by the user, or the value of ROOT returned on the
!    previous call to the function.
!
!    Input, real ( kind = rk ) FX, the value of the function at X.
!
!    Output, real ( kind = rk ) FERR, the smallest value of F encountered.
!
!    Output, real ( kind = rk ) XERR, the width of the change-in-sign interval,
!    if one was encountered.
!
!    Input/output, real ( kind = rk ) Q(9), storage needed by the function.
!    Before the first call, the user must set Q(1) to 0.
!
!    Output, real ( kind = rk ) ROOT_RC, an improved estimate for the root.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) d
  real ( kind = rk ) decr
  real ( kind = rk ) ferr
  real ( kind = rk ) fx
  integer i
  real ( kind = rk ) p
  real ( kind = rk ) q(9)
  real ( kind = rk ) r
  real ( kind = rk ) r8_huge
  real ( kind = rk ) r8_sign
  real ( kind = rk ) root_rc
  real ( kind = rk ) u
  real ( kind = rk ) v
  real ( kind = rk ) w
  real ( kind = rk ) x
  real ( kind = rk ) xerr
  real ( kind = rk ) z
!
!  If we found an exact zero, there is nothing more to do.
!
  if ( fx == 0.0D+00 ) then
    ferr = 0.0D+00
    xerr = 0.0D+00
    root_rc = x
    return
  end if

  ferr = abs ( fx )
!
!  If this is the first time, initialize, estimate the first root, and exit.
!
  if ( q(1) == 0.0D+00 ) then
    q(1) = fx
    q(2) = x
    q(3:9) = 0.0D+00
    root_rc = x + fx
    xerr = r8_huge ( )
    return
  end if
!
!  This is not the first call.
!
  q(9) = q(9) + 1.0D+00
!
!  Check for too many iterations.
!
  if ( 80.0D+00 < q(9) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROOT_RC - Fatal error!'
    write ( *, '(a,i6)' ) '  Number of iterations = ', int ( q(9) )
    stop
  end if
!
!  Check for a repeated X value.
!
  if ( ( 2.0 <= q(9) .and. x == q(4) ) .or. x == q(2) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROOT_RC - Fatal error!'
    write ( *, '(a,i6)' ) '  Value of X has been input before.'
    stop
  end if
!
!  Push X -> A -> B -> C
!
  do i = 6, 3, -1
    q(i) = q(i-2)
  end do
  q(1) = fx
  q(2) = x
!
!  If we have a change-in-sign interval, store the opposite value.
!
  if ( r8_sign ( q(1) ) /= r8_sign ( q(3) ) ) then
    q(7) = q(3)
    q(8) = q(4)
  end if
!
!  Calculate XERR.
!
  if ( q(7) /= 0.0D+00 ) then
    xerr = abs ( q(8) - q(2) )
  else
    xerr = r8_huge ( )
  end if
!
!  If more than 30 iterations, and we have change-in-sign interval, bisect.
!
  if ( 30.0D+00 < q(9) .and. q(7) /= 0.0D+00 ) then
    root_rc = q(2) + ( q(8) - q(2) ) / 2.0D+00
    return
  end if

  v = ( q(3) - q(1) ) / ( q(4) - q(2) )
!
!  If 3 or more points, try Muller.
!
  if ( q(5) /= 0.0D+00 ) then
    u = ( q(5) - q(3) ) / ( q(6) - q(4) )
    w = q(4) - q(2)
    z = ( q(6) - q(2) ) / w
    r = ( z + 1.0D+00 ) * v - u

    if ( r /= 0.0D+00 ) then
      p = 2.0D+00 * z * q(1) / r
      d = 2.0D+00 * p / ( w * r ) * ( v - u )
      if ( -1.0D+00 <= d ) then
        root_rc = q(2) - p / ( 1.0D+00 + sqrt ( 1.0D+00 + d ) )
        if ( q(7) == 0.0D+00 .or. &
             ( q(2) < root_rc .and. root_rc < q(8) ) .or. &
             ( q(8) < root_rc .and. root_rc < q(2) ) ) then
          return
        end if
      end if
    end if
  end if
!
!  Try the secant step.
!
  if ( q(1) /= q(3) .or. q(7) == 0.0D+00 ) then
    if ( q(1) == q(3) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ROOT_RC - Fatal error!'
      write ( *, '(a)' ) '  Cannot apply any method.'
      stop
    end if
    decr = q(1) / v
    if ( abs ( decr ) * 4.6D+18 < abs ( q(2) ) ) then
      decr = 1.74D-18 * abs ( q(2) ) * r8_sign ( decr )
    end if
    root_rc = q(2) - decr
    if ( q(7) == 0.0D+00 .or. &
        ( q(2) < root_rc .and. root_rc < q(8) ) .or. &
        ( q(8) < root_rc .and. root_rc < q(2) ) ) then
      return
    end if
  end if
!
!  Apply bisection.
!
  root_rc = q(2) + ( q(8) - q(2) ) / 2.0D+00

  return
end
function r8_epsilon ( )

!*****************************************************************************80
!
!! R8_EPSILON returns the R8 roundoff unit.
!
!  Discussion:
!
!    The roundoff unit is a number R which is a power of 2 with the
!    property that, to the precision of the computer's arithmetic,
!      1 < 1 + R
!    but
!      1 = ( 1 + R / 2 )
!
!    FORTRAN90 provides the superior library routine
!
!      EPSILON ( X )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) R8_EPSILON, the round-off unit.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_epsilon

  r8_epsilon = 2.220446049250313D-016

  return
end
function r8_huge ( )

!*****************************************************************************80
!
!! R8_HUGE returns a very large R8.
!
!  Discussion:
!
!    The value returned by this function is NOT required to be the
!    maximum representable R8.  This value varies from machine to machine,
!    from compiler to compiler, and may cause problems when being printed.
!    We simply want a "very large" but non-infinite number.
!
!    FORTRAN90 provides a built-in routine HUGE ( X ) that
!    can return the maximum representable number of the same datatype
!    as X, if that is what is really desired.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) R8_HUGE, a "huge" value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_huge

  r8_huge = 1.0D+30

  return
end
function r8_sign ( x )

!*****************************************************************************80
!
!! R8_SIGN returns the sign of an R8.
!
!  Discussion:
!
!    value = -1 if X < 0;
!    value = +1 if X => 0.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the number whose sign is desired.
!
!    Output, real ( kind = rk ) R8_SIGN, the sign of X:
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_sign
  real ( kind = rk ) x

  if ( x < 0.0D+00 ) then
    r8_sign = -1.0D+00
  else
    r8_sign = +1.0D+00
  end if

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
