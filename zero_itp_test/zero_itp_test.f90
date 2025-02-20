program main

!*****************************************************************************80
!
!! zero_itp_test() tests zero_itp().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 March 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) epsi
  real ( kind = rk ), external :: f_01
  real ( kind = rk ), external :: f_02
  real ( kind = rk ), external :: f_03
  real ( kind = rk ), external :: f_04
  real ( kind = rk ), external :: f_05
  real ( kind = rk ), external :: f_06
  real ( kind = rk ) k1
  real ( kind = rk ) k2
  integer n0
  character ( len = 80 ) title
  logical verbose

  write ( *, '(a)' ) ' '
  call timestamp ( )
  write ( *, '(a)' ) 'zero_itp_test():'
  write ( *, '(a)' ) '  Fortran90 version.'
  write ( *, '(a)' ) '  Test zero_itp(), which seeks a root of a function F(X)'
  write ( *, '(a)' ) '  in an interval [A,B].'

  a = 1.0D+00
  b = 2.0D+00
  epsi = sqrt ( epsilon ( 1.0D+00 ) )
  k1 = 1.0D+00 / ( b - a ) / 5.0D+00
  k2 = 2.0D+00
  n0 = 1
  verbose = .false.
  title = 'f_01(x) = sin ( x ) - x / 2'
  call zero_itp_example ( f_01, a, b, epsi, k1, k2, n0, verbose, title )

  a = 0.0D+00
  b = 1.0D+00
  epsi = sqrt ( epsilon ( 1.0D+00 ) )
  k1 = 1.0D+00 / ( b - a ) / 5.0D+00
  k2 = 2.0D+00
  n0 = 1
  verbose = .false.
  title = 'f_02(x) = 2 * x - exp ( - x )'
  call zero_itp_example ( f_02, a, b, epsi, k1, k2, n0, verbose, title )

  a = -1.0D+00
  b =  0.5D+00
  epsi = sqrt ( epsilon ( 1.0D+00 ) )
  k1 = 1.0D+00 / ( b - a ) / 5.0D+00
  k2 = 2.0D+00
  n0 = 1
  verbose = .false.
  title = 'f_03(x) = x * exp ( - x )'
  call zero_itp_example ( f_03, a, b, epsi, k1, k2, n0, verbose, title )

  a =  0.0001D+00
  b =  20.0D+00
  epsi = sqrt ( epsilon ( 1.0D+00 ) )
  k1 = 1.0D+00 / ( b - a ) / 5.0D+00
  k2 = 2.0D+00
  n0 = 1
  verbose = .false.
  title = 'f_04(x) = exp ( x ) - 1 / ( 100 * x * x )'
  call zero_itp_example ( f_04, a, b, epsi, k1, k2, n0, verbose, title )

  a = -5.0D+00
  b =  2.0D+00
  epsi = sqrt ( epsilon ( 1.0D+00 ) )
  k1 = 1.0D+00 / ( b - a ) / 5.0D+00
  k2 = 2.0D+00
  n0 = 1
  verbose = .false.
  title = 'f_05(x) = (x+3) * (x-1) * (x-1)'
  call zero_itp_example ( f_05, a, b, epsi, k1, k2, n0, verbose, title )

  a = 1.0D+00
  b =  2.0D+00
  epsi = 0.0005D+00
  k1 = 1.0D+00 / ( b - a ) / 5.0D+00
  k2 = 2.0D+00
  n0 = 1
  verbose = .false.
  title = 'f_06(x) = x^3 - x - 2'
  call zero_itp_example ( f_06, a, b, epsi, k1, k2, n0, verbose, title )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'zero_itp_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine zero_itp_example ( f, a, b, epsi, k1, k2, n0, verbose, title )

!*****************************************************************************80
!
!! zero_itp_example() tests zero_itp() on one test function.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 May 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    external real ( kind = rk ) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose zero is being sought.
!
!    real ( kind = rk ) A, B, the endpoints of the change of sign interval.
!
!    real ( kind = rk ) epsi: error tolerance between exact and computed roots.
!
!    real ( kind = rk ) k1: a parameter, with suggested value 0.2 / ( b - a ).
!
!    real ( kind = rk ) k2: a parameter, typically set to 2.
!
!    integer n0: a parameter that can be set to 0 for difficult problems,
!    but is usually set to 1, to take more advantage of the secant method.
!
!    logical verbose: if true, requests additional output from zero_itp().
!
!    character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer calls
  real ( kind = rk ) epsi
  real ( kind = rk ), external :: f
  real ( kind = rk ) fa
  real ( kind = rk ) fb
  real ( kind = rk ) fz
  real ( kind = rk ) k1
  real ( kind = rk ) k2
  integer n0
  character ( len = *  ) title
  logical verbose
  real ( kind = rk ) z
  
  call zero_itp ( f, a, b, epsi, k1, k2, n0, verbose, z, fz, calls )

  fa = f ( a )
  fb = f ( b )

  write ( *, '(a)' ) ' '
  write ( *, '(2x,a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      A                 Z             B'
  write ( *, '(a)' ) '    F(A)              F(Z)          F(B)'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,f14.8,2x,f14.8,2x,f14.8)' ) a,  z,  b
  write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) fa, fz, fb
  write ( *, '(a,i8)' ) '  Number of calls to F = ', calls
  write ( *, '(a,g14.6)' ) '  Tolerance epsi = ', epsi 
  write ( *, '(a,g14.6,a,g14.6,a,i4)' ) &
    '  Parameters k1 =', k1, ', k2 = ', k2, ', n0 = ', n0


  return
end
function f_01 ( x )

!*****************************************************************************80
!
!! f_01() evaluates sin ( x ) - x / 2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) X, the evaluation point.
!
!  Output:
!
!    real ( kind = rk ) F_01, the value of the function at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f_01
  real ( kind = rk ) x

  f_01 = sin ( x ) - 0.5D+00 * x

  return
end
function f_02 ( x )

!*****************************************************************************80
!
!! f_02() evaluates 2*x-exp(-x).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) X, the evaluation point.
!
!  Output:
!
!    real ( kind = rk ) F_02, the value of the function at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f_02
  real ( kind = rk ) x

  f_02 = 2.0D+00 * x - exp ( - x )

  return
end
function f_03 ( x )

!*****************************************************************************80
!
!! f_03() evaluates x*exp(-x).
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) X, the evaluation point.
!
!  Output:
!
!    real ( kind = rk ) F_03, the value of the function at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f_03
  real ( kind = rk ) x

  f_03 = x * exp ( - x )

  return
end
function f_04 ( x )

!*****************************************************************************80
!
!! f_04() evaluates exp(x) - 1 / (100*x*x).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) X, the evaluation point.
!
!  Output:
!
!    real ( kind = rk ) F_04, the value of the function at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f_04
  real ( kind = rk ) x

  f_04 = exp ( x ) - 1.0D+00 / 100.0D+00 / x / x

  return
end
function f_05 ( x )

!*****************************************************************************80
!
!! f_05() evaluates (x+3)*(x-1)*(x-1).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) X, the evaluation point.
!
!  Output:
!
!    real ( kind = rk ) F_05, the value of the function at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f_05
  real ( kind = rk ) x

  f_05 = ( x + 3.0D+00 ) * ( x - 1.0D+00 ) * ( x - 1.0D+00 )

  return
end
function f_06 ( x )

!*****************************************************************************80
!
!! f_06() evaluates x^3 - x - 2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 March 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) X, the evaluation point.
!
!  Output:
!
!    real ( kind = rk ) F_06, the value of the function at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f_06
  real ( kind = rk ) x

  f_06 = x ** 3 - x - 2.0D+00

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

