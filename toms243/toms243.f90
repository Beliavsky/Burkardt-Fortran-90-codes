subroutine toms243 ( z, value )

!*****************************************************************************80
!
!! toms243() computes the natural logarithm for complex values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 January 2018
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Collens,
!    Algorithm 243: Logarithm of a Complex Number,
!    Communications of the Association for Computing Machinery,
!    Volume 7, Number 11, November 1964, page 660.
!
!  Parameters:
!
!    Input, complex ( kind = ck ) Z, the argument of the function.
!
!    Output, complex ( kind = ck ) VALUE, the value of the function.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) e
  real ( kind = rk ) f
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  complex ( kind = ck ) value
  complex ( kind = ck ) z

  a = real ( z, kind = rk )
  b = aimag ( z )
!
!  Ugly hack to get NaN values.
!
  if ( a == 0.0D+00 .and. b == 0.0D+00 ) then
    c = 1.0D+00 / a
    d = 1.0D+00 / b
  else
    e = a / 2.0D+00
    f = b / 2.0D+00
    if ( abs ( e ) < 0.5D+00 .and. abs ( f ) < 0.5D+00 ) then
      c = abs ( 2.0D+00 * a ) + abs ( 2.0D+00 * b )
      d = 8.0D+00 * ( a / c ) * a + 8.0D+00 * ( b / c ) * b
      c = 0.5D+00 * ( log ( c ) + log ( d ) ) - log ( sqrt ( 8.0D+00 ) )
    else
      c = abs ( e / 2.0D+00 ) + abs ( f / 2.0D+00 )
      d = 0.5D+00 * ( e / c ) * e + 0.5D+00 * ( f / c ) * f
      c = 0.5D+00 * ( log ( c ) + log ( d ) ) + log ( sqrt ( 8.0D+00 ) )
    end if

    if ( ( a /= 0.0D+00 ) .and. abs ( f ) <= abs ( e ) ) then
      if ( sign ( 1.0D+00, a ) /= -1.0D+00 ) then
        d = atan ( b / a )
      else if ( sign ( 1.0D+00, b ) /= -1.0D+00 ) then
        d = atan ( b / a ) + r8_pi
      else
        d = atan ( b / a ) - r8_pi
      end if
    else
      d = - atan ( a / b ) + r8_pi / 2.0D+00 * sign ( 1.0D+00, b )
    end if

  end if

  value = cmplx ( c, d, kind = ck )

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
!  Parameters:
!
!    None
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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
