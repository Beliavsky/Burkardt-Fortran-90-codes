program main

!*****************************************************************************80
!
!! toms493_test() tests toms493().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'toms493_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test toms493().'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'toms493_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() uses a polynomial suggested by Driessen and Hunt.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 June 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    HB Driessen, EW Hunt,
!    Remark on Algorithm 429,
!    Communications of the ACM,
!    Volume 16, Number 9, page 579, September 1973.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: degree = 4
 
  real ( kind = rk8 ) c(degree+1)
  logical fail
  integer i
  real ( kind = rk8 ) pvi
  real ( kind = rk8 ) pvr
  real ( kind = rk8 ) root_i(degree)
  real ( kind = rk8 ) root_r(degree)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  p(x) = x^4 + 5.6562x^3 + 5.8854x^2'
  write ( *, '(a)' ) '             + 7.3646x   + 6.1354'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Approximate roots:'
  write ( *, '(a)' ) '    -1.0,'
  write ( *, '(a)' ) '    -4.674054017161709, '
  write ( *, '(a)' ) '     0.0089 + 1.1457 i,'
  write ( *, '(a)' ) '     0.0089 - 1.1457 i.'

  c(1) = 1.0D+00
  c(2) = 5.6562D+00
  c(3) = 5.8854D+00
  c(4) = 7.3646D+00
  c(5) = 6.1354D+00

  call rpoly ( c, degree, root_r, root_i, fail )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) &
    '      real(z)         imag(z)       real(f(z))      imag(f(z))'
  write ( *, '(a)' ) ''

  do i = 1, degree

    call polyev ( degree + 1, root_r(i), root_i(i), c, pvr, pvi )

    write ( *, '(2x,g18.10,2x,g18.10,2x,g14.6,2x,g14.6)' ) & 
      root_r(i), root_i(i), pvr, pvi

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! test02() tests a polynomial with a single repeated root.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: degree = 4

  real ( kind = rk8 ) c(degree+1)
  logical fail
  integer i
  real ( kind = rk8 ) pvi
  real ( kind = rk8 ) pvr
  real ( kind = rk8 ) root_i(degree)
  real ( kind = rk8 ) root_r(degree)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test02():'
  write ( *, '(a)' ) '  p(x) = x^4 - 8 x^3 + 24 x^2 - 32 x + 16'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Exact roots:'
  write ( *, '(a)' ) '    2 (multiplicity 4)'

  c(1) =   1.0D+00
  c(2) =  -8.0D+00
  c(3) =  24.0D+00
  c(4) = -32.0D+00
  c(5) =  16.0D+00

  call rpoly ( c, degree, root_r, root_i, fail )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) &
    '      real(z)         imag(z)       real(f(z))      imag(f(z))'
  write ( *, '(a)' ) ''

  do i = 1, degree

    call polyev ( degree + 1, root_r(i), root_i(i), c, pvr, pvi )

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      root_r(i), root_i(i), pvr, pvi

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! test03() tests a polynomial in the roots of unity.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: degree = 5

  real ( kind = rk8 ) c(degree+1)
  logical fail
  integer i
  real ( kind = rk8 ) pvi
  real ( kind = rk8 ) pvr
  real ( kind = rk8 ) root_i(degree)
  real ( kind = rk8 ) root_r(degree)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test03():'
  write ( *, '(a)' ) '  p(x) = x^5 - 1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Exact roots:'
  write ( *, '(a)' ) '    0.3090 + 0.9510 i'
  write ( *, '(a)' ) '   -0.8090 + 0.5877 i'
  write ( *, '(a)' ) '   -0.8090 - 0.5877 i'
  write ( *, '(a)' ) '    0.3090 - 0.9510 i'
  write ( *, '(a)' ) '    1'

  c(1) =  1.0D+00
  c(2) =  0.0D+00
  c(3) =  0.0D+00
  c(4) =  0.0D+00
  c(5) =  0.0D+00
  c(6) = -1.0D+00

  call rpoly ( c, degree, root_r, root_i, fail )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) &
    '      real(z)         imag(z)       real(f(z))      imag(f(z))'
  write ( *, '(a)' ) ''

  do i = 1, degree

    call polyev ( degree + 1, root_r(i), root_i(i), c, pvr, pvi )

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      root_r(i), root_i(i), pvr, pvi

  end do

  return
end
subroutine polyev ( nn, sr, si, c, pvr, pvi )

!*****************************************************************************80
!
!! polyev() evaluates a real polynomial by the Horner recurrence.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 June 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins, Joseph Traub,
!    Algorithm 419: Zeros of a Complex Polynomial,
!    Communications of the Association for Computing Machinery,
!    Volume 15, Number 2, February 1972, pages 97-99.
!
!  Input:
!
!    integer nn: the polynomial order.
!
!    real sr, si: the real and imaginary parts of the
!    polynomial argument.
!
!    real c(nn): the polynomial coefficients.
!
!  Output:
!
!    real pvr, pvi: the real and imaginary parts of the
!    polynomial value.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nn

  real ( kind = rk8 ) c(nn)
  integer i
  real ( kind = rk8 ) pvi
  real ( kind = rk8 ) pvr
  real ( kind = rk8 ) si
  real ( kind = rk8 ) sr
  real ( kind = rk8 ) t

  do i = 1, nn
    if ( i == 1 ) then
      pvr = c(1)
      pvi = 0.0
    else
      t =   pvr * sr - pvi * si + c(i)
      pvi = pvr * si + pvi * sr
      pvr = t
    end if
  end do

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
