subroutine llsq0 ( n, x, y, a )

!*****************************************************************************80
!
!! LLSQ solves a linear least squares problem matching y=a*x to data.
!
!  Discussion:
!
!    A formula for a line of the form Y = A * X is sought, which
!    will minimize the root-mean-square error to N data points ( X(I), Y(I) );
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 January 2019
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = rk ) X(N), Y(N), the coordinates of the data points.
!
!    Output, real ( kind = rk ) A the slope of the 
!    least-squares approximant to the data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ( kind = 4 ) n

  real ( kind = rk ) a
  real ( kind = rk ) bot
  real ( kind = rk ) top
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)
!
!  Special case.
!
  if ( n == 1 ) then
    if ( x(1) == 0.0D+00 ) then
      a = 1.0D+00
    else
      a = y(1) / x(1)
    end if
    return
  end if

  top = dot_product ( x(1:n), y(1:n) )
  bot = dot_product ( x(1:n), x(1:n) )

  a = top / bot

  return
end
subroutine llsq ( n, x, y, a, b )

!*****************************************************************************80
!
!! LLSQ solves a linear least squares problem matching a line to data.
!
!  Discussion:
!
!    A formula for a line of the form Y = A * X + B is sought, which
!    will minimize the root-mean-square error to N data points ( X(I), Y(I) );
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = rk ) X(N), Y(N), the coordinates of the data points.
!
!    Output, real ( kind = rk ) A, B, the slope and Y-intercept of the 
!    least-squares approximant to the data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ( kind = 4 ) n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) bot
  real ( kind = rk ) top
  real ( kind = rk ) x(n)
  real ( kind = rk ) xbar
  real ( kind = rk ) y(n)
  real ( kind = rk ) ybar
!
!  Special case.
!
  if ( n == 1 ) then
    a = 0.0D+00
    b = y(1)
    return
  end if
!
!  Average X and Y.
!
  xbar = sum ( x(1:n) ) / real ( n, kind = rk )
  ybar = sum ( y(1:n) ) / real ( n, kind = rk )
!
!  Compute Beta.
!
  top = dot_product ( x(1:n) - xbar, y(1:n) - ybar )
  bot = dot_product ( x(1:n) - xbar, x(1:n) - xbar )

  a = top / bot

  b = ybar - a * xbar

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
