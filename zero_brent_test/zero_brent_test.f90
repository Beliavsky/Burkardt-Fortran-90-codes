program main

!*****************************************************************************80
!
!! zero_brent_test() tests zero_brent().
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

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ), external :: f_01
  real ( kind = rk ), external :: f_02
  real ( kind = rk ), external :: f_03
  real ( kind = rk ), external :: f_04
  real ( kind = rk ), external :: f_05
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  call timestamp ( )
  write ( *, '(a)' ) 'zero_brent_test():'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test zero_brent(), which seeks a root of a function F(X)'
  write ( *, '(a)' ) '  in an interval [A,B].'

  a = 1.0D+00
  b = 2.0D+00
  title = 'f_01(x) = sin ( x ) - x / 2'
  call zero_example ( a, b, f_01, title )

  a = 0.0D+00
  b = 1.0D+00
  title = 'f_02(x) = 2 * x - exp ( - x )'
  call zero_example ( a, b, f_02, title )

  a = -1.0D+00
  b =  0.5D+00
  title = 'f_03(x) = x * exp ( - x )'
  call zero_example ( a, b, f_03 , title)

  a =  0.0001D+00
  b =  20.0D+00
  title = 'f_04(x) = exp ( x ) - 1 / ( 100 * x * x )'
  call zero_example ( a, b, f_04, title )

  a = -5.0D+00
  b =  2.0D+00
  title = 'f_05(x) = (x+3) * (x-1) * (x-1)'
  call zero_example ( a, b, f_05, title )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'zero_brent_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine zero_example ( a, b, f, title )

!*****************************************************************************80
!
!! zero_example() tests zero_brent() on one test function.
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
!    real ( kind = rk ) A, B, the endpoints of the change of sign interval.
!
!    external real ( kind = rk ) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose zero is being sought.
!
!    character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer calls
  real ( kind = rk ), external :: f
  real ( kind = rk ) fa
  real ( kind = rk ) fb
  real ( kind = rk ) fx
  real ( kind = rk ) t
  character ( len = *  ) title
  real ( kind = rk ) x
  real ( kind = rk ) zero_brent
  
  t = epsilon ( t )
  x = zero_brent ( a, b, t, f, calls )
  fx = f ( x )
  fa = f ( a )
  fb = f ( b )

  write ( *, '(a)' ) ' '
  write ( *, '(2x,a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      A                 X             B'
  write ( *, '(a)' ) '    F(A)              F(X)          F(B)'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,f14.8,2x,f14.8,2x,f14.8)' ) a,  x,  b
  write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) fa, fx, fb
  write ( *, '(a,i8)' ) '  Number of calls to F = ', calls

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

