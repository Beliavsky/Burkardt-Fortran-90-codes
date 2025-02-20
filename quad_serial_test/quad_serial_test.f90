program main

!*****************************************************************************80
!
!! quad_serial_test() tests quad_serial().
!
!  Discussion:
!
!    This code is a candidate for parallelization.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 October 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) error
  real ( kind = rk8 ) exact
  external f
  real ( kind = rk8 ) f
  integer n
  real ( kind = rk8 ) q
  real ( kind = rk8 ) quad_serial
  real ( kind = rk8 ) wtime
  real ( kind = rk8 ) wtime1
  real ( kind = rk8 ) wtime2

  a =  0.0D+00
  b = 10.0D+00
  n = 10000000
  exact = 0.49936338107645674464D+00

  call timestamp ( )
  write ( *, * ) ' '
  write ( *, * ) 'quad_serial_test():'
  write ( *, * ) '  Fortran90 version'
  write ( *, * ) '  Estimate the integral of f(x) from A to B.'
  write ( *, * ) '  f(x) = 50 / ( pi * ( 2500 * x * x + 1 ) ).'
  write ( *, * ) ' '
  write ( *, * ) '  A        = ', a
  write ( *, * ) '  B        = ', b
  write ( *, * ) '  N        = ', n

  call cpu_time ( wtime1 )
  q = quad_serial ( f, a, b, n )
  call cpu_time ( wtime2 )

  error = abs ( q - exact )
  wtime = wtime2 - wtime1
 
  write ( *, * ) ' '
  write ( *, * ) '  Exact    = ', exact
  write ( *, * ) '  Estimate = ', q
  write ( *, * ) '  Error    = ', error
  write ( *, * ) '  Time     = ', wtime
!
!  Terminate.
!
  write ( *, * ) ' '
  write ( *, * ) 'quad_serial_test():'
  write ( *, * ) '  Normal end of execution.'
  write ( *, * ) ' '
  call timestamp ( )

  stop
end
function f ( x )

!*****************************************************************************80
!
!! f() evaluates the function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk8 ) X, the evaluation point.
!
!  Output:
!
!    real ( kind = rk8 ) F, the function value.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) f
  real ( kind = rk8 ) pi
  real ( kind = rk8 ) x

  pi = 3.141592653589793D+00
  f = 50.0D+00 / ( pi * ( 2500.0D+00 * x * x + 1.0D+00 ) )

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

