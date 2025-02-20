program main

!*****************************************************************************80
!
!! orbital_decay() models orbital decay in a two body planetary system.
!
!  Discussion:
!
!    Given two massive bodies subject to gravity, it is possible to write down
!    differential equations describing their motion.  These equations are
!    simpler to formulate in the frame of reference in which the center of 
!    mass of the two bodies does not move.  If one body is much more massive 
!    than the other, then our calculations in this new frame are essentially
!    the same as in the original geometry.  This is the case when one body
!    is the sun, and another a planet.  
!
!    This simulation would need to be modified if we wanted to consider
!    the behavior of two bodies of comparable mass, and expected to see
!    them both moving, or, even in the sun-planet case, if we wanted to
!    allow the sun to have a velocity while we stayed in a fixed frame
!    of reference.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 May 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer neqn
  integer step_num
  real ( kind = rk8 ), allocatable :: ts(:)
  real ( kind = rk8 ), allocatable :: ys(:,:)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'orbital_decay():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  This simulation follows a small body for 20 orbits'
  write ( *, '(a)' ) '  around a relatively massive body - such as Mercury around'
  write ( *, '(a)' ) '  the sun.'
  write ( *, '(a)' ) '  Kepler''s equations for a two body system are used.'
  write ( *, '(a)' ) '  Initially, the orbit is NOT an ellipse, but as time passes,'
  write ( *, '(a)' ) '  the orbit decays into an elliptical shape.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use rkf45() for the ODE integrator.'

  neqn = 4
  step_num = 1000

  allocate ( ts(0:step_num) )
  allocate ( ys(neqn,0:step_num) )

  call rkf45_solve ( neqn, step_num, ts, ys )
!
!  Create graphics files for processing by gnuplot.
!
  call gnuplot_files ( neqn, step_num, ts, ys )
!
!  Free memory.
!
  deallocate ( ts )
  deallocate ( ys )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'orbital_decay()'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine gnuplot_files ( neqn, step_num, ts, ys )

!*****************************************************************************80
!
!! GNUPLOT_FILES creates two files for processing by gnuplot.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer NEQN, the number of equations.
!
!    integer STEP_NUM, the number of steps to take.
!
!    real ( kind = rk8 ) TS(0:STEP_NUM), the time values.
!
!    real ( kind = rk8 ) YS(NEQN,0:STEP_NUM), the solution values.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer neqn
  integer step_num 

  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data_filename
  integer data_unit
  integer j
  real ( kind = rk8 ) ts(0:step_num)
  real ( kind = rk8 ) ys(neqn,0:step_num)

  call get_unit ( data_unit )
  data_filename = 'orbital_decay_data.txt'
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do j = 0, step_num
    write ( data_unit, '(2x,i6,2x,5(2x,g14.6))' ) &
      j, ts(j), ys(1:neqn,j)
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) '  Created data file "' // trim ( data_filename ) // '".'

  call get_unit ( command_unit )
  command_filename = 'orbital_decay_commands.txt'
  open ( unit = command_unit, file = command_filename, status = 'replace' )
  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "orbital_decay.png"'
  write ( command_unit, '(a)' ) 'set xlabel "X"'
  write ( command_unit, '(a)' ) 'set ylabel "Y"'
  write ( command_unit, '(a)' ) 'set title "Orbital decay after 20 orbits"'
  write ( command_unit, '(a)' ) 'set size ratio -1'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'set style fill solid'
  write ( command_unit, '(a)' ) 'set object 1 circle fc rgb "red"'
  write ( command_unit, '(a)' ) 'set object 1 circle at 0,0 size 0.05'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 3:5 lw 3 linecolor rgb "blue"'

  write ( command_unit, '(a)' ) 'quit'
  close ( unit = command_unit )
  write ( *, '(a)' ) &
    '  Created command file "' // trim ( command_filename ) // '".'

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
subroutine rkf45_solve ( neqn, step_num, ts, ys )

!*****************************************************************************80
!
!! RKF45_SOLVE runs the two body ODE system.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 April 2011
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer NEQN, the number of equations.
!
!    integer STEP_NUM, the number of steps to take.
!
!  Output:
!
!    real ( kind = rk8 ) TS(0:STEP_NUM), the time values.
!
!    real ( kind = rk8 ) YS(NEQN,0:STEP_NUM), the solution values.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer neqn
  integer step_num 

  real ( kind = rk8 ) abserr
  integer flag
  external kepler
  real ( kind = rk8 ) relerr
  integer step
  real ( kind = rk8 ) t
  real ( kind = rk8 ) t_out
  real ( kind = rk8 ) t_start
  real ( kind = rk8 ) t_stop
  real ( kind = rk8 ) ts(0:step_num)
  real ( kind = rk8 ) y(neqn)
  real ( kind = rk8 ) yp(neqn)
  real ( kind = rk8 ) ys(neqn,0:step_num)

  abserr = 1.0D-10
  relerr = 1.0D-10

  flag = 1

  t_start = 0.0D+00
  t_stop = 20.0D+00 * 3.895D+00

  t = 0.0D+00
  t_out = 0.0D+00

  y(1:neqn) = (/ 1.0D+00, 0.0D+00,  0.0D+00,  0.50D+00 /)

  call kepler ( t, y, yp )

  ys(1:neqn,0) = y(1:neqn)
  ts(0) = t

  do step = 1, step_num

    t = ( real ( step_num - step + 1, kind = rk8 ) * t_start &
        + real (            step - 1, kind = rk8 ) * t_stop ) &
        / real ( step_num,            kind = rk8 )

    t_out = ( real ( step_num - step, kind = rk8 ) * t_start &
            + real (            step, kind = rk8 ) * t_stop ) &
            / real ( step_num,        kind = rk8 )

    call rkf45 ( kepler, neqn, y, yp, t, t_out, relerr, abserr, flag )

    if ( abs ( flag ) /= 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RKF45_SOLVE - Warning!'
      write ( *, '(a,i4,a,g14.6)' ) '  Output value of FLAG = ', flag, &
        ' at T_OUT = ', t_out
    end if

    ys(1:neqn,step) = y(1:neqn)
    ts(step) = t_out

  end do

  return
end
subroutine kepler ( t, u, up )

!*****************************************************************************80
!
!! KEPLER evaluates the right hand side of the Kepler ODE system.
!
!  Discussion:
!
!    The Kepler ODE system has the form
!
!      u' = kepler ( t, u )
!
!    where u is a vector of length 4 whose components are the position
!    and velocity of a small body orbiting a massive one.
!
!      u = [ x(t), x'(t), y(t), y'(t) ]
!      
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 May 2013
!
!  Input:
!
!    real ( kind = rk8 ) T, the current time.
!
!    real ( kind = rk8 ) U(4), the current state.
!
!  Output:
!
!    real ( kind = rk8 ) UP(4), the derivative of the current state.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) r3
  real ( kind = rk8 ) t
  real ( kind = rk8 ) u(4)
  real ( kind = rk8 ) up(4)

  call r8_fake_use ( t )

  r3 = sqrt ( ( u(1) ** 2 + u(3) ** 2 ) ** 3 )

  up = (/ u(2), -u(1) / r3, u(4), -u(3) / r3 /)

  return
end
subroutine r8_fake_use ( x )

!*****************************************************************************80
!
!! r8_fake_use() pretends to use an R8 variable.
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
!    15 August 2021
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

