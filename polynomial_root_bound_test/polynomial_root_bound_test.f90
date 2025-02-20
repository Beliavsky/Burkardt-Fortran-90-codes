program main

!*****************************************************************************80
!
!! polynomial_root_bound_test() tests polynomial_root_bound().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 December 2023
!
!  Author:
!
!    John Burkardt
!
  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'polynomial_root_bound_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test polynomial_root_bound()'

  call polynomial1_test ( )
  call polynomial2_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'polynomial_root_bound_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine polynomial1_test ( )

!*****************************************************************************80
!
!! polynomial1_test() deals with a particular polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 December 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) b
  complex ( kind = ck ) c(0:n)
  integer i
  real ( kind = rk ) polynomial_root_bound
  complex ( kind = ck ) r(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'polynomial1_test():'
  write ( *, '(a)' ) '  Bound the roots of a specific polynomial:'
  write ( *, '(a)' ) '    12z^5 + 2z^2 + 23i.'

  c = (/ (0.0,23.0), (0.0,0.0), (2.0,0.0), (0.0,0.0), (0.0,0.0), (12.0,0.0) /)
  
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Polynomial coefficients are:'
  do i = 0, n
    write ( *, '(2x,i2,2x,g14.6,2x,g14.6)' ) i, c(i)
  end do

  b = polynomial_root_bound ( n, c )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Root magnitude bound is ', b

  r = (/ &
    (-1.10401598, -0.33686293), (-0.66152059, 0.89701442), &
    ( 0.02570877, -1.13896251), ( 0.67740946, 0.94586564), &
    ( 1.06241834, -0.36705461) /)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Polynomial roots and norms are:'
  do i = 1, n
    write ( *, '(2x,i2,2x,g14.6,2x,g14.6,2x,g14.6)' ) i, r(i), abs ( r(i) )
  end do

  return
end
subroutine polynomial2_test ( )

!*****************************************************************************80
!
!! polynomial2_test() deals with a random polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 December 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 5

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) b
  complex ( kind = ck ) c(0:n)
  integer i
  real ( kind = rk ) polynomial_root_bound
  complex ( kind = ck ) r(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'polynomial2_test():'
  write ( *, '(a)' ) '  Bound the roots of a random polynomial:'

  call c8vec_normal_01 ( n, r )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Polynomial roots and norms are:'
  do i = 1, n
    write ( *, '(2x,i2,2x,g14.6,2x,g14.6,2x,g14.6)' ) i, r(i), abs ( r(i) )
  end do

  call roots_to_c8poly ( n, r, c )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Polynomial coefficients are:'
  do i = 0, n
    write ( *, '(2x,i2,2x,g14.6,2x,g14.6)' ) i, c(i)
  end do

  b = polynomial_root_bound ( n, c )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Root magnitude bound is ', b

  return
end
subroutine c8vec_normal_01 ( n, x )

!*****************************************************************************80
!
!! c8vec_normal_01() returns a unit pseudonormal C8VEC.
!
!  Discussion:
!
!    A C8VEC is an array of double precision complex values.
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 December 2023
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of values desired.
!
!  Output:
!
!    complex ( kind = ck ) X(N), a sample of the standard normal PDF.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) r(n)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  complex ( kind = ck ) x(n)

  call random_number ( harvest = r(1:n) )

  x(1:n) = sqrt ( - 2.0D+00 * log ( r(1:n) ) ) * &
    cmplx ( cos ( 2.0D+00 * r8_pi * r(1:n) ), &
            sin ( 2.0D+00 * r8_pi * r(1:n) ), kind = ck )

  return
end
subroutine roots_to_c8poly ( n, r, c )

!*****************************************************************************80
!
!! roots_to_c8poly() converts polynomial roots to polynomial coefficients.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of roots specified.
!
!    complex ( kind = ck ) R(N), the roots.
!
!  Output:
!
!    complex ( kind = ck ) C(0:N), the coefficients of the polynomial.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer n

  complex ( kind = ck ) c(0:n)
  integer i
  integer j
  complex ( kind = ck ) r(n)
!
!  Initialize C to (0, 0, ..., 0, 1).
!  Essentially, we are setting up a divided difference table.
!
  c(0:n-1) = 0.0D+00
  c(n) = 1.0D+00
!
!  Convert to standard polynomial form by shifting the abscissas
!  of the divided difference table to 0.
!
  do j = 1, n
    do i = 1, n + 1 - j
      c(n-i) = c(n-i) - r(n+1-i-j+1) * c(n-i+1)
    end do
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

