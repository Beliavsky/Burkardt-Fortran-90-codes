program main

!*****************************************************************************80
!
!! midpoint_fixed_test() tests midpoint_fixed().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 April 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  real ( kind = rk ) tspan(2)
  real ( kind = rk ) y0(2)

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'midpoint_fixed_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test midpoint_fixed() on several ODE''s.'

  tspan(1) = 0.0D+00
  tspan(2) = 5.0D+00
  y0(1) = 5000.0D+00
  y0(2) = 100.0D+00
  n = 200
  call predator_prey_midpoint_fixed_test ( tspan, y0, n )

  tspan(1) = 0.0D+00
  tspan(2) = 1.0D+00
  y0(1) = 0.0D+00
  n = 27
  call stiff_midpoint_fixed_test ( tspan, y0(1), n )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'midpoint_fixed_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine predator_prey_deriv ( t, rf, drfdt )

!*****************************************************************************80
!
!! predator_prey_deriv() evaluates the right hand side of the system.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 February 2020
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    George Lindfield, John Penny,
!    Numerical Methods Using MATLAB,
!    Second Edition,
!    Prentice Hall, 1999,
!    ISBN: 0-13-012641-1,
!    LC: QA297.P45.
!
!  Input:
!
!    real ( kind = rk ) T, the current time.
!
!    real ( kind = rk ) RF(2), the current solution variables, rabbits and foxes.
!
!  Output:
!
!    real ( kind = rk ) DRFDT(2), the right hand side of the 2 ODE's.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) drfdt(2)
  real ( kind = rk ) rf(2)
  real ( kind = rk ) t

  call r8_fake_use ( t )

  drfdt(1) =    2.0 * rf(1) - 0.001 * rf(1) * rf(2)
  drfdt(2) = - 10.0 * rf(2) + 0.002 * rf(1) * rf(2)

  return
end
subroutine predator_prey_midpoint_fixed_test ( tspan, p0, n )

!*****************************************************************************80
!
!! predator_prey_midpoint_fixed_test(): predator-prey using midpoint_fixed().
!
!  Discussion:
!
!    The physical system under consideration is a pair of animal populations.
!
!    The PREY reproduce rapidly for each animal alive at the beginning of the
!    year, two more will be born by the end of the year.  The prey do not have
!    a natural death rate instead, they only die by being eaten by the predator.
!    Every prey animal has 1 chance in 1000 of being eaten in a given year by
!    a given predator.
!
!    The PREDATORS only die of starvation, but this happens very quickly.
!    If unfed, a predator will tend to starve in about 1/10 of a year.
!    On the other hand, the predator reproduction rate is dependent on
!    eating prey, and the chances of this depend on the number of available prey.
!
!    The resulting differential equations can be written:
!
!      PREY(0) = 5000         
!      PRED(0) =  100
!
!      d PREY / dT =    2 * PREY(T) - 0.001 * PREY(T) * PRED(T)
!      d PRED / dT = - 10 * PRED(T) + 0.002 * PREY(T) * PRED(T)
!
!    Here, the initial values (5000,100) are a somewhat arbitrary starting point.
!
!    The pair of ordinary differential equations that result have an interesting
!    behavior.  For certain choices of the interaction coefficients (such as
!    those given here), the populations of predator and prey will tend to 
!    a periodic oscillation.  The two populations will be out of phase the number
!    of prey will rise, then after a delay, the predators will rise as the prey
!    begins to fall, causing the predator population to crash again.
!
!    There is a conserved quantity, which here would be:
!      E(r,f) = 0.002 r + 0.001 f - 10 ln(r) - 2 ln(f)
! 
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 February 2020
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    George Lindfield, John Penny,
!    Numerical Methods Using MATLAB,
!    Second Edition,
!    Prentice Hall, 1999,
!    ISBN: 0-13-012641-1,
!    LC: QA297.P45.
!
!  Input:
!
!    real ( kind = rk ) TSPAN = [ T0, TMAX ], the initial and final times.
!    A reasonable value might be [ 0, 5 ].
!
!    real ( kind = rk ) P0 = [ PREY, PRED ], the initial number of prey and
!    predators.  A reasonable value might be [ 5000, 100 ].
!
!    integer N: the number of time steps.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 2
  integer n

  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data_filename
  integer data_unit
  character ( len = * ), parameter :: header = 'predator_prey_midpoint_fixed'
  integer i
  real ( kind = rk ) p0(m)
  real ( kind = rk ) pout(n+1,m)
  external predator_prey_deriv
  real ( kind = rk ) t(n+1)
  real ( kind = rk ) tspan(2)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'predator_prey_midpoint_fixed_test()'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  A pair of ordinary differential equations for a population'
  write ( *, '(a)' ) '  of predators and prey are solved using midpoint_fixed().'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  The exact solution shows periodic behavior, with a fixed'
  write ( *, '(a)' ) '  period and amplitude.'

  call midpoint_fixed ( predator_prey_deriv, tspan, p0, n, m, t, pout )
!
!  Create the data file.
!
  call get_unit ( data_unit )

  data_filename = header // '_data.txt'

  open ( unit = data_unit, file = data_filename, status = 'replace' )

  do i = 1, n + 1
    write ( data_unit, '(5(2x,g14.6))' ) t(i), pout(i,1), pout(i,2)
  end do

  close ( unit = data_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  predator_prey_midpoint_fixed_test: data stored in "' &
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
  write ( command_unit, '(a)' ) 'set output "' // trim ( header ) // '.png"'
  write ( command_unit, '(a)' ) 'set xlabel "<-- Prey -->"'
  write ( command_unit, '(a)' ) 'set ylabel "<-- Predator -->"'
  write ( command_unit, '(a)' ) 'set title "midpoint fixed: predator prey ODE"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 2:3 with lines lw 3'
  write ( command_unit, '(a)' ) 'quit'

  close ( unit = command_unit )

  write ( *, '(a)' ) '  predator_prey_midpoint_fixed_test: plot commands stored in "' &
    // trim ( command_filename ) // '".'

  return
end
subroutine stiff_deriv ( t, y, dydt )

!*****************************************************************************80
!
!! stiff_deriv() evaluates the right hand side of the stiff ODE.
!
!  Discussion:
!
!    y' = 50 * ( cos(t) - y )
!    y(0) = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) T, Y: the time and solution value.
!
!  Output:
!
!    real ( kind = rk ) DYDT: the derivative value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) dydt
  real ( kind = rk ) t
  real ( kind = rk ) y

  dydt = 50.0D+00 * ( cos ( t ) - y )

  return
end
subroutine stiff_exact ( n, t, y )

!*****************************************************************************80
!
!! stiff_exact() evaluates the exact solution of the stiff ODE.
!
!  Discussion:
!
!    y' = 50 * ( cos(t) - y )
!    y(0) = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N: the number of values.
!
!    real ( kind = rk ) T(N): the evaluation times.
!
!  Output:
!
!    real ( kind = rk ) Y(N): the exact solution values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) t(n)
  real ( kind = rk ) y(n)

  y(1:n) = 50.0D+00 * ( sin ( t ) + 50.0D+00 * cos(t) &
    - 50.0D+00 * exp ( - 50.0D+00 * t ) ) / 2501.0D+00

  return
end
subroutine stiff_midpoint_fixed_test ( tspan, y0, n )

!*****************************************************************************80
!
!! stiff_midpoint() uses the midpoint method on the stiff ODE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) TSPAN(2): the first and last times.
!
!    real ( kind = rk ) Y0: the initial condition.
!
!    integer N: the number of steps to take.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 1
  integer n
  integer, parameter :: n2 = 101

  external stiff_deriv
  real ( kind = rk ) t1(n+1)
  real ( kind = rk ) t2(n2)
  real ( kind = rk ) tspan(2)
  real ( kind = rk ) y1(n+1)
  real ( kind = rk ) y0
  real ( kind = rk ) y2(n2)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'stiff_midpoint_fixed_test'
  write ( *, '(a)' ) '  Solve stiff ODE using midpoint_fixed()'

  call midpoint_fixed ( stiff_deriv, tspan, y0, n, m, t1, y1 )

  call r8vec_linspace ( n2, tspan(1), tspan(2), t2 )
  call stiff_exact ( n2, t2, y2 )

  call plot2 ( n+1, t1, y1, n2, t2, y2, 'stiff_midpoint_fixed', &
    'Stiff ODE: midpoint fixed' )

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! get_unit() returns a free FORTRAN unit number.
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
!  Output:
!
!    integer IUNIT, the free unit number.
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
subroutine plot2 ( n1, t1, y1, n2, t2, y2, header, title )

!*****************************************************************************80
!
!! plot2() plots two curves together.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N1: the size of the first data set.
!
!    real ( kind = rk ) T1(N1), Y1(N1), the first dataset.
!
!    integer N2: the size of the second data set.
!
!    real ( kind = rk ) T2(N2), Y2(N2), the secod dataset.
!
!    character ( len = * ) HEADER: an identifier for the data.
!
!    character ( len = * ) TITLE: a title to appear in the plot.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n1
  integer n2

  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data1_filename
  character ( len = 255 ) data2_filename
  integer data_unit
  character ( len = * ) header
  integer i
  character ( len = * ) title
  real ( kind = rk ) t1(n1)
  real ( kind = rk ) t2(n2)
  real ( kind = rk ) y1(n1)
  real ( kind = rk ) y2(n2)
!
!  Create the data files.
!
  call get_unit ( data_unit )
  data1_filename = header // '_data1.txt'
  open ( unit = data_unit, file = data1_filename, status = 'replace' )
  do i = 1, n1
    write ( data_unit, '(5(2x,g14.6))' ) t1(i), y1(i)
  end do
  close ( unit = data_unit )

  call get_unit ( data_unit )
  data2_filename = header // '_data2.txt'
  open ( unit = data_unit, file = data2_filename, status = 'replace' )
  do i = 1, n2
    write ( data_unit, '(5(2x,g14.6))' ) t2(i), y2(i)
  end do
  close ( unit = data_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  plot2: data stored in "' &
    // trim ( data1_filename ) // '" and "' // trim ( data2_filename ) // '".'
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
  write ( command_unit, '(a)' ) 'set output "' // trim ( header ) // '.png"'
  write ( command_unit, '(a)' ) 'set xlabel "<-- T -->"'
  write ( command_unit, '(a)' ) 'set ylabel "<-- Y(T) -->"'
  write ( command_unit, '(a)' ) 'set title "' // trim ( title ) // '"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data1_filename ) // &
    '" using 1:2 with lines lw 3 lt rgb "red",\'
  write ( command_unit, '(a)' ) '     "' // trim ( data2_filename ) // &
    '" using 1:2 with lines lw 3 lt rgb "blue"'
  write ( command_unit, '(a)' ) 'quit'

  close ( unit = command_unit )

  write ( *, '(a)' ) '  plot2: plot commands stored in "' &
    // trim ( command_filename ) // '".'

  return
end
subroutine r8_fake_use ( x )

!*****************************************************************************80
!
!! r8_fake_use() pretends to use a variable.
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
!    real ( kind = rk ) X, the variable to be "used".
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) x

  if ( x /= x ) then
    write ( *, '(a)' ) '  r8_fake_use: variable is NAN.'
  end if

  return
end

subroutine r8vec_linspace ( n, a, b, x )

!*****************************************************************************80
!
!! r8vec_linspace() creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
!
!    In other words, the interval is divided into N-1 even subintervals,
!    and the endpoints of intervals are used as the points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of entries in the vector.
!
!    real ( kind = rk ) A, B, the first and last entries.
!
!  Output:
!
!    real ( kind = rk ) X(N), a vector of linearly spaced data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer i
  real ( kind = rk ) x(n)

  if ( n == 1 ) then

    x(1) = ( a + b ) / 2.0D+00

  else

    do i = 1, n
      x(i) = ( real ( n - i,     kind = rk ) * a   &
             + real (     i - 1, kind = rk ) * b ) &
             / real ( n     - 1, kind = rk )
    end do

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
 
