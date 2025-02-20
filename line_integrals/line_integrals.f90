function line01_length ( )

!*****************************************************************************80
!
!! line01_length(): length of the unit line in 1D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) LINE01_LENGTH, the length.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) line01_length

  line01_length = 1.0D+00

  return
end
subroutine line01_monomial_integral ( e, integral )

!*****************************************************************************80
!
!! LINE01_MONOMIAL_INTEGRAL: monomial integral over the unit line in 1D.
!
!  Discussion:
!
!    The integration region is 
!
!      0 <= X <= 1.
!
!    The monomial is F(X) = X^E.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Academic Press, 1984, page 263.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) E, the exponent.
!    E must not equal -1.
!
!    Output, real ( kind = rk ) INTEGRAL, the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ( kind = 4 ) e
  real ( kind = rk ) integral

  if ( e == -1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE01_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  Exponent E = -1 is not allowed!'
    stop 1
  end if

  integral = 1.0D+00 / real ( e + 1, kind = rk )

  return
end
subroutine line01_sample ( n, x )

!*****************************************************************************80
!
!! LINE01_SAMPLE samples the unit line in 1D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Output, real ( kind = rk ) X(N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ( kind = 4 ) n

  real ( kind = rk ) x(n)

  call random_number ( harvest = x(1:n) )

  return
end
subroutine monomial_value_1d ( n, e, x, value )

!*****************************************************************************80
!
!! MONOMIAL_VALUE_1D evaluates a monomial in 1D.
!
!  Discussion:
!
!    This routine evaluates a monomial of the form
!
!      x^expon
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) E, the exponent.
!
!    Input, real ( kind = rk ) X(N), the point coordinates.
!
!    Output, real ( kind = rk ) VALUE(N), the value of the monomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ( kind = 4 ) n

  integer ( kind = 4 ) e
  real ( kind = rk ) value(n)
  real ( kind = rk ) x(n)

  value(1:n) = x(1:n) ** e

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

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

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
