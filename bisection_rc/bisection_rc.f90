subroutine bisection_rc ( a, b, x, fx, job )

!*****************************************************************************80
!
!! bisection_rc() seeks a zero of f(x) in a change of sign interval.
!
!  Discussion:
!
!    The bisection method is used.
!
!    This routine uses reverse communication, so that the function is always
!    evaluated in the calling program.
!
!    On the first call, the user sets JOB = 0, and the values of A and B.
!    Thereafter, the user checks the returned value of JOB and follows 
!    directions.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) A, B, the endpoints of the change of sign interval.
!
!    real ( kind = rk ) FX, the function value at X.
!
!    integer JOB, a communication flag.
!    The user sets JOB to 0 before the first call.
!
!  Output:
!
!    real ( kind = rk ) A, B, the updated endpoints of the change of sign interval.
!
!    real ( kind = rk ) X, a point at which the function is to be evaluated.
!
!    integer JOB: indicates the action required of the user:
!    1: please evaluate F(X) at the output value X.  If |F(X)| is not
!    small enough, call this function again.
!  
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ), save :: fa
  real ( kind = rk ), save :: fb
  real ( kind = rk ) fx
  integer job
  real ( kind = rk ) r8_sign
  integer, save :: state
  real ( kind = rk ) x

  if ( job == 0 ) then

    fa = 0.0D+00
    fb = 0.0D+00
    state = 1
    x = a
    job = 1

  else if ( state == 1 ) then

    fa = fx
    x = b
    state = 2

  else if ( state == 2 ) then

    fb = fx

    if ( r8_sign ( fa ) == r8_sign ( fb ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BISECTION_RC - Fatal error!'
      write ( *, '(a)' ) '  F(A) and F(B) have the same sign.'
      stop
    end if

    x = ( a + b ) / 2.0D+00
    state = 3

  else

    if ( r8_sign ( fx ) == r8_sign ( fa ) ) then
      a = x
      fa = fx
    else
      b = x
      fb = fx
    end if
    x = ( a + b ) / 2.0D+00
    state = 3

  end if

  return
end
function r8_sign ( x )

!*****************************************************************************80
!
!! r8_sign() returns the sign of an R8.
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
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) X, the number whose sign is desired.
!
!  Output:
!
!    real ( kind = rk ) R8_SIGN, the sign of X:
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
!! timestamp() prints the current YMDHMS date as a time stamp.
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
!    03 September 2021
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
