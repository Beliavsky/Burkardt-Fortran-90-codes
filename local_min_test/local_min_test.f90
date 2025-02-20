program main

!*****************************************************************************80
!
!! local_min_test() tests local_min().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 June 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ), external :: g_01
  real ( kind = rk ), external :: g_02
  real ( kind = rk ), external :: g_03
  real ( kind = rk ), external :: g_04
  real ( kind = rk ), external :: g_05
  real ( kind = rk ), external :: g_06
  real ( kind = rk ), external :: g_07

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'local_min_test():'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test local_min(), which seeks'
  write ( *, '(a)' ) '  a local minimizer of a function F(X)'
  write ( *, '(a)' ) '  in an interval [A,B].'

  a = 0.0D+00
  b = 3.141592653589793D+00
  call local_min_example ( a, b, g_01, &
    'g_01(x) = ( x - 2 ) * ( x - 2 ) + 1' )

  a = 0.0D+00
  b = 1.0D+00
  call local_min_example ( a, b, g_02, &
    'g_02(x) = x * x + exp ( - x )' )

  a = -2.0D+00
  b =  2.0D+00
  call local_min_example ( a, b, g_03, &
    'g_03(x) = x^4 + 2x^2 + x + 3' )

  a =  0.0001D+00
  b =  1.0D+00
  call local_min_example ( a, b, g_04, &
    'g_04(x) = exp ( x ) + 1 / ( 100 x )' )

  a =  0.0002D+00
  b = 2.0D+00
  call local_min_example ( a, b, g_05, &
    'g_05(x) = exp ( x ) - 2x + 1/(100x) - 1/(1000000x^2)' )

  a = 1.8D+00
  b = 1.9D+00
  call local_min_example ( a, b, g_06, &
    'g_06(x) = -x*sin(10*pi*x)-1.0' )

  a = -1.2D+00
  b =  2.7D+00
  call local_min_example ( a, b, g_07, &
    'g_07(x) = max(-2(x-1),8(x-1)) + 25*(x-1)^2' )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'local_min_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine local_min_example ( a, b, f, title )

!*****************************************************************************80
!
!! local_min_example() tests local_min() on one test function.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 May 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) A, B, the endpoints of the interval.
!
!    external real ( kind = rk ) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose local minimum is being sought.
!
!    character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer ( kind = 4 ) calls
  real ( kind = rk ) eps
  real ( kind = rk ), external :: f
  real ( kind = rk ) fa
  real ( kind = rk ) fb
  real ( kind = rk ) fx
  real ( kind = rk ) local_min
  real ( kind = rk ) t
  character ( len = * ) title
  real ( kind = rk ) x

  eps = epsilon ( eps )
  t = sqrt ( eps )

  fx = local_min ( a, b, eps, t, f, x, calls )
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
function g_01 ( x )

!*****************************************************************************80
!
!! g_01() evaluates (x-2)^2 + 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 April 2008
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
!    real ( kind = rk ) G_01, the value of the function at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) g_01
  real ( kind = rk ) x

  g_01 = ( x - 2.0D+00 ) * ( x - 2.0D+00 ) + 1.0D+00

  return
end
function g_02 ( x )

!*****************************************************************************80
!
!! g_02() evaluates x^2 + exp ( - x ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 April 2008
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
!    real ( kind = rk ) G_02, the value of the function at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) g_02
  real ( kind = rk ) x

  g_02 = x * x + exp ( - x )

  return
end
function g_03 ( x )

!*****************************************************************************80
!
!! g_03() evaluates x^4+2x^2+x+3.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 April 2008
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
!    real ( kind = rk ) G_03, the value of the function at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) g_03
  real ( kind = rk ) x

  g_03 = ( ( x * x + 2.0D+00 ) * x + 1.0D+00 ) * x + 3.0D+00

  return
end
function g_04 ( x )

!*****************************************************************************80
!
!! g_04() evaluates exp(x)+1/(100X)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 April 2008
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
!    real ( kind = rk ) G_04, the value of the function at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) g_04
  real ( kind = rk ) x

  g_04 = exp ( x ) + 0.01D+00 / x

  return
end
function g_05 ( x )

!*****************************************************************************80
!
!! g_05() evaluates exp(x) - 2x + 1/(100x) - 1/(1000000x^2)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 April 2008
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
!    real ( kind = rk ) G_05, the value of the function at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) g_05
  real ( kind = rk ) x

  g_05 = exp ( x ) - 2.0D+00 * x + 0.01D+00 / x - 0.000001D+00 / x / x

  return
end
function g_06 ( x )

!*****************************************************************************80
!
!! g_06() evaluates - x * sin(10 pi x ) - 1.0;
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 May 2021
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
!    real ( kind = rk ) G_06, the value of the function at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) g_06
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) x

  g_06 = - x * sin ( 10.0D+00 * r8_pi * x ) - 1.0D+00

  return
end
function g_07 ( x )

!*****************************************************************************80
!
!! g_07() evaluates max(-2(x-1), 8(x-1)) + 25 (x-1)^2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 May 2021
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
!    real ( kind = rk ) G_07, the value of the function at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) g_07
  real ( kind = rk ) x

  g_07 = max ( -2.0D+00 * ( x - 1.0D+00 ), &
                8.0D+00 * ( x - 1.0D+00 ) ) &
       + 25.0D+00 * ( x - 1.0D+00 ) ** 2

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

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

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

