program main

!*****************************************************************************80
!
!! logistic_exact_test() tests logistic_exact().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 May 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  interface
    subroutine logistic_parameters ( r_in, k_in, t0_in, y0_in, tstop_in, &
      r_out, k_out, t0_out, y0_out, tstop_out )
      integer, parameter :: rk8 = kind ( 1.0D+00 )
      real ( kind = rk8 ), optional :: r_in
      real ( kind = rk8 ), optional :: r_out
      real ( kind = rk8 ), optional :: k_in
      real ( kind = rk8 ), optional :: k_out
      real ( kind = rk8 ), optional :: t0_in
      real ( kind = rk8 ), optional :: t0_out
      real ( kind = rk8 ), optional :: y0_in
      real ( kind = rk8 ), optional :: y0_out
      real ( kind = rk8 ), optional :: tstop_in
      real ( kind = rk8 ), optional :: tstop_out
    end subroutine
  end interface

  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data_filename
  integer data_unit
  character ( len = 14 ) header
  integer i
  integer j
  real ( kind = rk8 ) k
  integer n
  real ( kind = rk8 ) r
  real ( kind = rk8 ), allocatable :: t(:)
  real ( kind = rk8 ) t0
  real ( kind = rk8 ) tstop
  real ( kind = rk8 ), allocatable :: y(:,:)
  real ( kind = rk8 ) y0
  real ( kind = rk8 ), save, dimension ( 4 ) :: y0_vec = (/ &
    0.05D+00, 0.1D+00, 0.2D+00, 0.5D+00 /)

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'logistic_exact_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  logistic_exact() evaluates exact solutions of'
  write ( *, '(a)' ) '  the logistic ODE.'

  header = 'logistic_exact'
!
!  Get parameters.
!
  call logistic_parameters ( r_out = r, k_out = k, t0_out = t0, y0_out = y0, &
    tstop_out = tstop )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  parameters:'
  write ( *, '(a,g14.6)' ) '    r     = ', r
  write ( *, '(a,g14.6)' ) '    k     = ', k
  write ( *, '(a,g14.6)' ) '    t0    = ', t0
  write ( *, '(a,g14.6)' ) '    y0    = ', y0
  write ( *, '(a,g14.6)' ) '    tstop = ', tstop

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Evaluate exact solution for varied values of y0:'

  n = 101
  allocate ( t(1:n) )
  allocate ( y(1:n,1:4) )

  call r8vec_linspace ( n, t0, tstop, t )

  do j = 1, 4

    y0 = y0_vec(j)
    write( *, '(a,g14.6)' ) '    y0 = ', y0
    call logistic_parameters ( y0_in = y0 )
    call logistic_exact ( n, t, y(1:n,j) )

  end do
!
!  Create the data file.
!
  call get_unit ( data_unit )

  data_filename = header // '_data.txt'

  open ( unit = data_unit, file = data_filename, status = 'replace' )

  do i = 1, n
    write ( data_unit, '(5(2x,g14.6))' ) t(i), y(i,1:4)
  end do

  close ( unit = data_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data stored in "' // trim ( data_filename ) // '".'
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
  write ( command_unit, '(a)' ) 'set xlabel "<-- t -->"'
  write ( command_unit, '(a)' ) 'set ylabel "<-- y(t) -->"'
  write ( command_unit, '(a)' ) 'set title "flame_exact:"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 1:2  with lines lw 3,\'
  write ( command_unit, '(a)' ) '     "' // trim ( data_filename ) // &
    '" using 1:3 with lines lw 3,\'
  write ( command_unit, '(a)' ) '     "' // trim ( data_filename ) // &
    '" using 1:4 with lines lw 3,\'
  write ( command_unit, '(a)' ) '     "' // trim ( data_filename ) // &
    '" using 1:5 with lines lw 3'
  write ( command_unit, '(a)' ) 'quit'

  close ( unit = command_unit )

  write ( *, '(a)' ) '  Plot commands stored in "' // trim ( command_filename ) // '".'
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'logistic_exact_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

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
!    free FORTRAN unit.  The code assumes that units 5 and 6
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
