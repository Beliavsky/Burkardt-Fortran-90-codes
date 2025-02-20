program main

!*****************************************************************************80
!
!! sigmoid_derivative_test() tests sigmoid_derivative().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 May 2024
!
!  Author:
!
!    John Burkardt
!
  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'sigmoid_derivative_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test sigmoid_derivative.'

  call sigmoid_derivative_coef_test ( )
  call sigmoid_derivative_value_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'sigmoid_derivative_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine sigmoid_derivative_coef_test ( )

!*****************************************************************************80
!
!! sigmoid_derivative_coef_test() tests sigmoid_derivative_coef().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 May 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ), allocatable :: coef(:)
  character ( len = 80 ) label
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'sigmoid_derivative_coef_test():'
  write ( *, '(a)' ) '  sigmoid_derivative_coef() returns the coefficients of'
  write ( *, '(a)' ) '  the expansion of the nth derivative of the sigmoid'
  write ( *, '(a)' ) '  function in terms of powers of the sigmoid function.'

  do n = 0, 4

    allocate ( coef(1:n+2) )
    call sigmoid_derivative_coef ( n, coef )
    write ( label, '(a,i1,a)' ) '  s^(', n, ')(x)'
    call sigmoid_poly_print ( n + 2, coef, label )
    deallocate ( coef )

  end do

  return
end
subroutine sigmoid_derivative_value_test ( )

!*****************************************************************************80
!
!! sigmoid_derivative_value_test() tests sigmoid_derivative_value().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 May 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: nvec = 51

  character ( len = 255 ) header
  integer i
  integer n
  real ( kind = rk8 ) sigmoid_derivative_value
  real ( kind = rk8 ) xvec(nvec)
  real ( kind = rk8 ) yvec(nvec)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'sigmoid_derivative_value_test():'
  write ( *, '(a)' ) '  sigmoid_derivative_value() evaluates the nth derivative'
  write ( *, '(a)' ) '  of the sigmoid function at the location x.'

  do n = 0, 3

    call r8vec_linspace ( nvec, -5.0D+00, +5.0D+00, xvec )
    yvec(1:nvec) = 0.0D+00
    do i = 1, nvec
      yvec(i) = sigmoid_derivative_value ( n, xvec(i) )
    end do

    write ( header, '(a,i1)' ) 'sigmoid_derivative_', n
    call sigmoid_derivative_plot ( nvec, xvec, yvec, trim ( header ) )

  end do

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
subroutine sigmoid_derivative_plot ( nvec, xvec, yvec, header )

!*****************************************************************************80
!
!! sigmoid_derivative_plot() plots a sigmoid derivative.
!
!  Discussion:
!
!    Actually, we simply create two files for processing by gnuplot().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 May 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nvec

  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data_filename
  integer data_unit
  character ( len = * ) header
  integer i
  character ( len = 255 ) png_filename
  real ( kind = rk8 ) xvec(nvec)
  real ( kind = rk8 ) yvec(nvec)
!
!  Create a graphics data file.
!
  data_filename = header // '_data.txt'
  call get_unit ( data_unit )

  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do i = 1, nvec
    write ( data_unit, '(2g14.6)' ) xvec(i), yvec(i)
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Created graphics data file "' &
    // trim ( data_filename ) // '".'

  png_filename = header // '.png'
!
!  Create graphics command file.
!
  call get_unit ( command_unit )
  command_filename = header // '_commands.txt'
  open ( unit = command_unit, file = command_filename, status = 'replace' )
  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "' // trim ( png_filename ) // '"'
  write ( command_unit, '(a)' ) 'set xlabel "<--- X --->"'
  write ( command_unit, '(a)' ) 'set ylabel "<-- Y(X) -->"'
  write ( command_unit, '(a)' ) 'set title "' // header // '"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 1:2 lw 3 linecolor rgb "blue",\'
  close ( unit = command_unit )
  write ( *, '(a)' ) &
    '  Created command file "' // trim ( command_filename ) // '".'

  return
end
subroutine sigmoid_poly_print ( n, a, title )

!*****************************************************************************80
!
!! sigmoid_poly_print() prints a polynomial in s(x).
!
!  Discussion:
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1) * x + ... + a(n-1) * x^(n-1) + a(n) * x^(n)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the dimension of A.
!
!    real ( kind = rk ) A(0:N), the polynomial coefficients.
!    A(0) is the constant term and
!    A(N) is the coefficient of X^N.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(0:n-1)
  integer i
  real ( kind = rk ) mag
  character plus_minus
  character ( len = * ) title

  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) trim ( title ) // ' = '

  do i = 0, n - 1

    if ( a(i) < 0.0D+00 ) then
      plus_minus = '-'
    else
      plus_minus = '+'
    end if

    mag = abs ( a(i) )

    if ( mag /= 0.0D+00 ) then

      if ( 2 <= i ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * s(x) ^ '', i3 )' ) &
          plus_minus, mag, i
      else if ( i == 1 ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * s(x)'' )' ) plus_minus, mag
      else if ( i == 0 ) then
        write ( *, ' ( ''         '', a1, g14.6 )' ) plus_minus, mag
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
!    real ( kind = rk8 ) A, B, the first and last entries.
!
!  Output:
!
!    real ( kind = rk8 ) X(N), a vector of linearly spaced data.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  integer i
  real ( kind = rk8 ) x(n)

  if ( n == 1 ) then

    x(1) = ( a + b ) / 2.0D+00

  else

    do i = 1, n
      x(i) = ( real ( n - i,     kind = rk8 ) * a   &
             + real (     i - 1, kind = rk8 ) * b ) &
             / real ( n     - 1, kind = rk8 )
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
