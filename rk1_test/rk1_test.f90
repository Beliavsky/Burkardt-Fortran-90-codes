program main

!*****************************************************************************80
!
!! rk1_test() tests rk1().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 November 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer n

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'rk1_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test euler().'

  n = 50
  call rk1_humps_test ( n )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'rk1_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop ( 0 )
end
subroutine rk1_humps_test ( n )

!*****************************************************************************80
!
!! rk1_humps_test() calls rk1() for the humps problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 November 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N: the number of steps to take.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: m = 1
  integer n

  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data_filename
  integer data_unit
  character ( len = * ), parameter :: header = 'rk1_humps'
  external humps_deriv
  integer i
  real ( kind = rk8 ) t(n+1)
  real ( kind = rk8 ) tspan(2)
  real ( kind = rk8 ) y(n+1,1)
  real ( kind = rk8 ) y0(m)
  real ( kind = rk8 ) y2(n+1,1)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'rk1_humps_test():'

  tspan(1) = 0.0D+00
  tspan(2) = 2.0D+00
  call humps_exact ( 1, tspan(1), y0 )

  call rk1 ( humps_deriv, tspan, y0, n, m, t, y )

  call humps_exact ( n+1, t, y2 )
!
!  Create the data file.
!
  call get_unit ( data_unit )

  data_filename = header // '_data.txt'

  open ( unit = data_unit, file = data_filename, status = 'replace' )

  do i = 1, n + 1
    write ( data_unit, '(3(2x,g14.6))' ) t(i), y(i,1), y2(i,1)
  end do

  close ( unit = data_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  rk1_humps_test(): data stored in "' &
    // trim ( data_filename ) // '".'
!
!  Create the command file.
!
  call get_unit ( command_unit )

  command_filename = trim ( header ) // '_commands.txt'

  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "' // trim ( header ) // '_test.png"'
  write ( command_unit, '(a)' ) 'set xlabel "<-- T -->"'
  write ( command_unit, '(a)' ) 'set ylabel "<-- Y(T) -->"'
  write ( command_unit, '(a)' ) 'set title "rk1: Humps ODE"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 1:2 with lines lw 3 lt rgb "red", \'
  write ( command_unit, '(a)' ) '     "' // trim ( data_filename ) // &
    '" using 1:3 with lines lw 3 lt rgb "blue"'

  write ( command_unit, '(a)' ) 'quit'

  close ( unit = command_unit )

  write ( *, '(a)' ) '  rk1_humps_test(): plot commands stored in "' &
    // trim ( command_filename ) // '".'

  return
end
subroutine humps_deriv ( x, y, yp )

!*****************************************************************************80
!
!! humps_deriv evaluates the right hand side of the humps ODE.
!
!  Discussion:
!
!    y = 1.0 / ( ( x - 0.3 )^2 + 0.01 )
!      + 1.0 / ( ( x - 0.9 )^2 + 0.04 )
!      - 6.0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 November 2023
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk8 ) x, y(1): the argument.
!
!  Output:
!
!    real ( kind = rk8 ) yp(1): the value of the derivative at x.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) x
  real ( kind = rk8 ) y(1)
  real ( kind = rk8 ) yp(1)

  call r8_fake_use ( y(1) )

  yp(1) = - 2.0D+00 * ( x - 0.3D+00 ) / ( ( x - 0.3D+00 )**2 + 0.01D+00 )**2 &
          - 2.0D+00 * ( x - 0.9D+00 ) / ( ( x - 0.9D+00 )**2 + 0.04D+00 )**2

  return
end
subroutine humps_exact ( n, x, y )

!*****************************************************************************80
!
!! humps_exact evaluates the humps function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n: the number of evaluation points.
!
!    real ( kind = rk8 ) x(n): the evaluation points.
!
!  Output:
!
!    real ( kind = rk8 ) y(n): the function values.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) x(n)
  real ( kind = rk8 ) y(n)

  write ( *, * ) 'x(1) = ', x(1)
  y = 1.0D+00 / ( ( x - 0.3D+00 )**2 + 0.01D+00 ) &
    + 1.0D+00 / ( ( x - 0.9D+00 )**2 + 0.04D+00 ) &
    - 6.0D+00

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer i
  integer ios
  integer iunit
  logical lopen
 
  iunit = 0
 
  do i = 1, 99
 
    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )
 
      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if
 
  end do

  return
end
subroutine r8_fake_use ( x )

!*****************************************************************************80
!
!! r8_fake_use pretends to use a variable.
!
!  Discussion:
!
!    Some compilers will issue a warning if a variable is unused.
!    Sometimes there's a good reason to include a variable in a program,
!    but not to use it.  Calling this function with that variable as
!    the argument will shut the compiler up.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk8 ) X, the variable to be "used".
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) x

  if ( x /= x ) then
    write ( *, '(a)' ) '  r8_fake_use(): variable is NAN.'
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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

