program main

!*****************************************************************************80
!
!! toms419_test() tests toms419().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'toms419_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test toms419().'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
!
!  Termination.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'toms419_test():'
  write ( *, '(a)' ) '  Normal end of execution.'

  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() tests the first polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: degree = 10
  integer, parameter :: order = degree + 1

  logical fail
  real ( kind = rk8 ) fzr(degree)
  real ( kind = rk8 ) fzi(degree)
  integer i
  real ( kind = rk8 ) pr(order)
  real ( kind = rk8 ) pi(order)
  real ( kind = rk8 ) zr(degree)
  real ( kind = rk8 ) zi(degree)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  Polynomial with zeros 1,2,...,10.'
  write ( *, '(a,i2)' ) '  Degree ', degree
 
  pr(1) = 1.0
  pr(2) = -55.0
  pr(3) = 1320.0
  pr(4) = -18150.0
  pr(5) = 157773.0
  pr(6) = -902055.0
  pr(7) = 3416930.0
  pr(8) = -8409500.0
  pr(9) = 12753576.0
  pr(10) = -10628640.0
  pr(11) = 3628800.0
  do i = 1, order
    pi(i) = 0.0
  end do

  call print_coefficients ( order, pr, pi )
  call cpoly ( pr, pi, degree, zr, zi, fail )

  if ( fail ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  cpoly() failed on this example.'
  else
    call print_zeros ( degree, zr, zi )
    call compute_fz ( degree, zr, zi, pr, pi, fzr, fzi )
    call print_fz ( degree, fzr, fzi )
  end if

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! test02() tests the second polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: degree = 3
  integer, parameter :: order = degree + 1

  logical fail
  real ( kind = rk8 ) fzr(degree)
  real ( kind = rk8 ) fzi(degree)
  real ( kind = rk8 ) pr(order)
  real ( kind = rk8 ) pi(order)
  real ( kind = rk8 ) zr(degree)
  real ( kind = rk8 ) zi(degree)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test02():'
  write ( *, '(a)' ) '  Zeros on imaginary axis.'
  write ( *, '(a,i2)' ) '  Degree ', degree

  pr(1) = 1.0
  pr(2) = 0.0
  pr(3) = -10001.0001D0
  pr(4) = 0.0
  pi(1) = 0.0
  pi(2) = -10001.0001d0
  pi(3) = 0.0
  pi(4) = 1.0

  call print_coefficients ( order, pr, pi )
  call cpoly ( pr, pi, degree, zr, zi, fail )

  if ( fail ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  cpoly() failed on this example.'
  else
    call print_zeros ( degree, zr, zi )
    call compute_fz ( degree, zr, zi, pr, pi, fzr, fzi )
    call print_fz ( degree, fzr, fzi )
  end if

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! test03() tests the third polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: degree = 10
  integer, parameter :: order = degree + 1

  logical fail
  real ( kind = rk8 ) fzr(degree)
  real ( kind = rk8 ) fzi(degree)
  real ( kind = rk8 ) pr(order)
  real ( kind = rk8 ) pi(order)
  real ( kind = rk8 ) zr(degree)
  real ( kind = rk8 ) zi(degree)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test03():'
  write ( *, '(a)' ) '  zeros at 1+I,1/2*(1+I)....1/(2**-9)*(1+I)'
  write ( *, '(a,i2)' ) '  Degree ', degree

  pr(1) = 1.0
  pr(2) = -1.998046875
  pr(3) = 0.0
  pr(4) = 0.7567065954208374D0
  pr(5) = -0.2002119533717632D0
  pr(6) = 1.271507365163416D-2
  pr(7) = 0.0
  pr(8) = -1.154642632172909D-5
  pr(9) = 1.584803612786345D-7
  pr(10) = -4.652065399568528D-10
  pr(11) = 0.0
  pi(1) = 0.0
  pi(2) = pr(2)
  pi(3) = 2.658859252929688d0
  pi(4) = -7.567065954208374d-1
  pi(5) = 0.0
  pi(6) = pr(6)
  pi(7) = -7.820779428584501d-4
  pi(8) = -pr(8)
  pi(9) = 0.0
  pi(10) = pr(10)
  pi(11) = 9.094947017729282d-13

  call print_coefficients ( order, pr, pi )
  call cpoly ( pr, pi, degree, zr, zi, fail )

  if ( fail ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  cpoly() failed on this example.'
  else
    call print_zeros ( degree, zr, zi )
    call compute_fz ( degree, zr, zi, pr, pi, fzr, fzi )
    call print_fz ( degree, fzr, fzi )
  end if

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! test04() tests the fourth polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: degree = 10
  integer, parameter :: order = degree + 1

  logical fail
  real ( kind = rk8 ) fzr(degree)
  real ( kind = rk8 ) fzi(degree)
  real ( kind = rk8 ) pr(order)
  real ( kind = rk8 ) pi(order)
  real ( kind = rk8 ) zr(degree)
  real ( kind = rk8 ) zi(degree)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test04():'
  write ( *, '(a)' ) '  Multiple zeros.'
  write ( *, '(a,i2)' ) '  Degree ', degree

  pr(1) = 1.0
  pr(2) = -10.0
  pr(3) = 3.0
  pr(4) = 284.0
  pr(5) = -1293.0
  pr(6) = 2374.0
  pr(7) = -1587.0
  pr(8) = -920.0
  pr(9) = 2204.0
  pr(10) = -1344.0
  pr(11) = 288.0
  pi(1) = 0.0
  pi(2) = -10.0
  pi(3) = 100.0
  pi(4) = -334.0
  pi(5) = 200.0
  pi(6) = 1394.0
  pi(7)  = -3836.0
  pi(8) = 4334.0
  pi(9) = -2352.0
  pi(10) = 504.0
  pi(11) = 0.0

  call print_coefficients ( order, pr, pi )
  call cpoly ( pr, pi, degree, zr, zi, fail )

  if ( fail ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  cpoly() failed on this example.'
  else
    call print_zeros ( degree, zr, zi )
    call compute_fz ( degree, zr, zi, pr, pi, fzr, fzi )
    call print_fz ( degree, fzr, fzi )
  end if

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! test05() tests the fifth polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: degree = 12
  integer, parameter :: order = degree + 1

  logical fail
  real ( kind = rk8 ) fzr(degree)
  real ( kind = rk8 ) fzi(degree)
  real ( kind = rk8 ) pr(order)
  real ( kind = rk8 ) pi(order)
  real ( kind = rk8 ) zr(degree)
  real ( kind = rk8 ) zi(degree)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test05():'
  write ( *, '(a)' ) '  Zeros on unit circle centered at 0+2i.'
  write ( *, '(a,i2)' ) '  Degree ', degree

  pr(1) = 1.0
  pr(2) = 0.0
  pr(3) = -264.0
  pr(4) = 0.0
  pr(5) = 7920.0
  pr(6) = 0.0
  pr(7) = -59136.0
  pr(8) = 0.0
  pr(9) = 126720.0
  pr(10) = 0.0
  pr(11) = -67584.0
  pr(12) = 0.0
  pr(13) = 4095.0
  pi(1) = 0.0
  pi(2) = -24.0
  pi(3) = 0.0
  pi(4) = 1760.0
  pi(5) = 0.0
  pi(6) = -25344.0
  pi(7) = 0.0
  pi(8) = 101376.0
  pi(9) = 0.0
  pi(10) = -112640.0
  pi(11) = 0.0
  pi(12) = 24576.0
  pi(13) = 0.0

  call print_coefficients ( order, pr, pi )
  call cpoly ( pr, pi, degree, zr, zi, fail )

  if ( fail ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  cpoly() failed on this example.'
  else
    call print_zeros ( degree, zr, zi )
    call compute_fz ( degree, zr, zi, pr, pi, fzr, fzi )
    call print_fz ( degree, fzr, fzi )
  end if

  return
end
subroutine print_coefficients ( order, pr, pi )

!*****************************************************************************80
!
!! print_coefficients() prints the coefficients.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
!  Input:
!
!    integer order: the order of the polynomial.
!
!    real ( kind = rk8 ) pr(order), pi(order): the real and imaginary parts
!    of the coefficients.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer order

  integer i
  real ( kind = rk8 ) pi(order)
  real ( kind = rk8 ) pr(order)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Coefficients:'
  write ( *, '(50(2d26.16/))' ) ( pr(i), pi(i), i = 1, order )

  return
end
subroutine print_zeros ( degree, zr, zi )

!*****************************************************************************80
!
!! print_zeros() prints the zeros.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
!  Input:
!
!    integer degree: the degree of the polynomial.
!
!    real ( kind = rk8 ) zr(degree), zi(degree): the real and imaginary parts
!    of the zeros.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer degree

  integer i
  real ( kind = rk8 ) zi(degree)
  real ( kind = rk8 ) zr(degree)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Zeros:'
  write ( *, '(50(2d26.16/))' ) ( zr(i), zi(i), i = 1, degree )

  return
end
subroutine compute_fz ( degree, zr, zi, pr, pi, fzr, fzi )

!*****************************************************************************80
!
!! compute_fz() computes the function values at the zeros.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
!  Input:
!
!    integer degree: the degree of the polynomial.
!
!    real ( kind = rk8 ) pr(degree+1), pi(degree+1): the real and imaginary 
!    parts of the polynomial coefficients.
!
!    real ( kind = rk8 ) zr(degree), zi(degree): the real and imaginary parts
!    of the computed zeros.
!
!  Output:
!
!    real ( kind = rk8 ) fzr(degree), fzi(degree): the real and imaginary parts
!    of the function values.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer degree

  integer i
  real ( kind = rk8 ) fzi(degree)
  real ( kind = rk8 ) fzr(degree)
  real ( kind = rk8 ) pi(degree+1)
  real ( kind = rk8 ) pr(degree+1)
  real ( kind = rk8 ) qpi(degree+1)
  real ( kind = rk8 ) qpr(degree+1)
  real ( kind = rk8 ) zi(degree+1)
  real ( kind = rk8 ) zr(degree+1)

  do i = 1, degree
    call polyev ( degree + 1, zr(i), zi(i), pr, pi, qpr, qpi, fzr(i), fzi(i) )
  end do

  return
end
subroutine print_fz ( degree, fzr, fzi )

!*****************************************************************************80
!
!! print_fz() prints the function values at the zeros.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
!  Input:
!
!    integer degree: the degree of the polynomial.
!
!    real ( kind = rk8 ) fzr(degree), fzi(degree): the real and imaginary parts
!    of the function values.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer degree

  integer i
  real ( kind = rk8 ) fzi(degree)
  real ( kind = rk8 ) fzr(degree)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Function values:'
  write ( *, '(50(2d26.16/))' ) ( fzr(i), fzi(i), i = 1, degree )

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
