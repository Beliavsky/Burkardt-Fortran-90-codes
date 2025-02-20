program main

!*****************************************************************************80
!
!! zero_laguerre_test() tests zero_laguerre().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 March 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'zero_laguerre_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test zero_laguerre().'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'zero_laguerre_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() runs the tests on a polynomial function.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) abserr
  integer degree
  real ( kind = rk ), external :: func01
  integer ierror
  integer k
  integer kmax
  real ( kind = rk ) x
  real ( kind = rk ) x0
!
!  Give a starting point.
!
  x0 = 1.0D+00
!
!  The polynomial degree.
!
  degree = 3
!
!  Set the error tolerance.
!
  abserr = 0.00001D+00
!
!  KMAX is the maximum number of iterations.
!
  kmax = 30

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  (Polynomial function F(X))'
  write ( *, '(a)' ) '  Find a root of F(X)=(X+3)*(X+3)*(X-2)=0'

  call laguerre ( x0, degree, abserr, kmax, func01, x, ierror, k )

  write ( *, '(a)' ) ' '
  if ( ierror /= 0 ) then
    write ( *, '(a,i2)' ) '  Iteration failed with ierror = ', ierror
  else
    write ( *, '(a,i2)' ) '  Iteration steps taken: ', k
    write ( *, '(a,g14.6)' ) '  Estimated root X = ', x
    write ( *, '(a,g14.6)' ) '  F(X) = ', func01 ( x, 0 )
  end if

  return
end
function func01 ( x, ider )

!*****************************************************************************80
!
!! func01() computes the function value for the first test.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real X, the point at which the evaluation is to take place.
!
!    integer IDER, specifies what is to be evaluated:
!    0, evaluate the function.
!    1, evaluate the first derivative.
!    2, evaluate the second derivative.
!    3, evaluate the third derivative.
!
!  Output:
!
!    real FUNC01, the value of the function or derivative.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) func01
  integer ider
  real ( kind = rk ) x

  if ( ider == 0 ) then
    func01 = ( x + 3.0D+00 )**2 * ( x - 2.0E+00 )
  else if ( ider == 1 ) then
    func01 = ( x + 3.0D+00 ) * ( 3.0E+00 * x - 1.0E+00 )
  else if ( ider == 2 ) then
    func01 = 6.0D+00 * x + 8.0E+00
  else if ( ider == 3 ) then
    func01 = 6.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'func01(): Fatal error!'
    write ( *, '(a,i8)' ) '  Derivative of order IDER = ', ider
    write ( *, '(a)' ) '  was requested.'
    stop
  end if

  return
end




subroutine test02 ( )

!*****************************************************************************80
!
!! test02() runs the tests on the Newton polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 March 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) abserr
  integer degree
  real ( kind = rk ), external :: func02
  integer ierror
  integer k
  integer kmax
  real ( kind = rk ) x
  real ( kind = rk ) x0
!
!  Give a starting point.
!
  x0 = 1.0D+00
!
!  The polynomial degree.
!
  degree = 3
!
!  Set the error tolerance.
!
  abserr = 0.00001D+00
!
!  KMAX is the maximum number of iterations.
!
  kmax = 30

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test02():'
  write ( *, '(a)' ) '  (Polynomial function F(X))'
  write ( *, '(a)' ) '  p(x) = x^3 - 2x - 5'

  call laguerre ( x0, degree, abserr, kmax, func02, x, ierror, k )

  write ( *, '(a)' ) ' '
  if ( ierror /= 0 ) then
    write ( *, '(a,i2)' ) '  Iteration failed with ierror = ', ierror
  else
    write ( *, '(a,i2)' ) '  Iteration steps taken: ', k
    write ( *, '(a,g14.6)' ) '  Estimated root X = ', x
    write ( *, '(a,g14.6)' ) '  F(X) = ', func02 ( x, 0 )
  end if

  return
end
function func02 ( x, ider )

!*****************************************************************************80
!
!! func02() computes the function value for the Newton polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 March 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real X, the point at which the evaluation is to take place.
!
!    integer IDER, specifies what is to be evaluated:
!    0, evaluate the function.
!    1, evaluate the first derivative.
!    2, evaluate the second derivative.
!    3, evaluate the third derivative.
!
!  Output:
!
!    real FUNC02, the value of the function or derivative.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) func02
  integer ider
  real ( kind = rk ) x

  if ( ider == 0 ) then
    func02 = x**3 - 2.0D+00 * x - 5.0D+00
  else if ( ider == 1 ) then
    func02 = 3.0D+00 * x**2 - 2.0D+00
  else if ( ider == 2 ) then
    func02 = 6.0D+00 * x 
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'func02(): Fatal error!'
    write ( *, '(a,i8)' ) '  Derivative of order IDER = ', ider
    write ( *, '(a)' ) '  was requested.'
    stop
  end if

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

