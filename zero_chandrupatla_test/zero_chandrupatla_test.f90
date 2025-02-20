program main

!*****************************************************************************80
!
!! zero_chandrupatla_test() tests zero_chandrupatla().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 March 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ), external :: f_01
  real ( kind = rk8 ), external :: f_02
  real ( kind = rk8 ), external :: f_03
  real ( kind = rk8 ), external :: f_04
  real ( kind = rk8 ), external :: f_05
  real ( kind = rk8 ), external :: f_06
  real ( kind = rk8 ), external :: f_07
  real ( kind = rk8 ), external :: f_08
  real ( kind = rk8 ), external :: f_09
  character ( len = 100 ) title

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'zero_chandrupatla_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  zero_chandrupatla() seeks a root of a function f(x)'
  write ( *, '(a)' ) '  in an interval [a,b].'

  a = 2.0
  b = 3.0
  title = 'f_01(x) = x**3 - 2 x - 5'
  call zero_chandrupatla_example ( f_01, a, b, title )

  a = 0.5
  b = 1.51
  title = 'f_02(x) = 1 - 1/x**2'
  call zero_chandrupatla_example ( f_02, a, b, title )

  a = 0.0
  b = 5.0
  title = 'f_03(x) = ( x - 3 )**3'
  call zero_chandrupatla_example ( f_03, a, b, title )

  a = 0.0
  b = 5.0
  title = 'f_04(x) = 6 * ( x - 2 )**5'
  call zero_chandrupatla_example ( f_04, a, b, title )

  a = -1.0
  b = 4.0
  title = 'f_05(x) = x**9'
  call zero_chandrupatla_example ( f_05, a, b, title )

  a = -1.0
  b = 4.0
  title = 'f_06(x) = x**19'
  call zero_chandrupatla_example ( f_06, a, b, title )

  a = -1.0
  b = 4.0
  title = 'f_07(x) = x e**(-1/x2)'
  call zero_chandrupatla_example ( f_07, a, b, title )

  a = 0.0002
  b = 2.0
  title = 'f_08(x) = -(3062(1-xi)e**(-x)/(xi+(1-xi)e**(-x)) - 1013 + 1628/x'
  call zero_chandrupatla_example ( f_08, a, b, title )

  a = 0.0002
  b = 1.0
  title = 'f_09(x) = e**x - 2 - 0.01/x**2 + 0.000002/x**3'
  call zero_chandrupatla_example ( f_09, a, b, title )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'zero_chandrupatla_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine zero_chandrupatla_example ( f, a, b, title )

!*****************************************************************************80
!
!! zero_chandrupatla_example() tests zero_chandrupatla() on a test function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 March 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    function f ( x ): the user-supplied function.
!
!    real a, b: the endpoints of the change of sign interval.
!
!    string title: a title for the problem.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  integer calls
  real ( kind = rk8 ), external :: f
  real ( kind = rk8 ) fa
  real ( kind = rk8 ) fb
  real ( kind = rk8 ) fz
  character ( len = * ) title
  real ( kind = rk8 ) z

  call zero_chandrupatla ( f, a, b, z, fz, calls )

  fz = f ( z )
  fa = f ( a )
  fb = f ( b )

  write ( *, '(a)' ) ''
  write ( *, '(2x,a)' ) trim ( title )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '      A                 Z             B'
  write ( *, '(a)' ) '    F(A)              F(Z)          F(B)'
  write ( *, '(a)' ) ''
  write ( *, '(2x,f14.6,2x,f14.6,2x,f14.6)' ) a,  z,  b
  write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) fa, fz, fb
  write ( *, '(a,i4)' ) '  Number of calls to F = ', calls

  return
end
function f_01 ( x )

!*****************************************************************************80
!
!! f_01() evaluates function 1.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) f_01
  real ( kind = rk8 ) x

  f_01 = x**3 - 2.0 * x - 5.0

  return
end
function f_02 ( x )

!*****************************************************************************80
!
!! f_02() evaluates function 2.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) f_02
  real ( kind = rk8 ) x

  f_02 = 1.0 - 1.0 / x**2

  return
end
function f_03 ( x )

!*****************************************************************************80
!
!! f_03() evaluates function 3.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) f_03
  real ( kind = rk8 ) x

  f_03 = ( x - 3.0 )**3

  return
end
function f_04 ( x )

!*****************************************************************************80
!
!! f_04() evaluates function 4.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) f_04
  real ( kind = rk8 ) x

  f_04 = 6.0 * ( x - 2.0 )**5

  return
end
function f_05 ( x )

!*****************************************************************************80
!
!! f_05() evaluates function 5.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) f_05
  real ( kind = rk8 ) x

  f_05 = x**9

  return
end
function f_06 ( x )

!*****************************************************************************80
!
!! f_06() evaluates function 6.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) f_06
  real ( kind = rk8 ) x

  f_06 = x**19

  return
end
function f_07 ( x )

!*****************************************************************************80
!
!! f_07() evaluates function 7.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) f_07
  real ( kind = rk8 ) x

  if ( abs ( x ) < 3.8D-04 ) then
    f_07 = 0.0
  else
    f_07 = x * exp ( - ( 1.0 / x**2 ) )
  end if

  return
end
function f_08 ( x )

!*****************************************************************************80
!
!! f_08() evaluates function 8.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) bot
  real ( kind = rk8 ) f_08
  real ( kind = rk8 ) top
  real ( kind = rk8 ) x
  real ( kind = rk8 ) xi

  xi = 0.61489
  top = ( 3062.0 * ( 1.0 - xi ) * exp ( - x ) )
  bot = ( xi + ( 1.0 - xi ) * exp ( - x ) )

  f_08 = - top / bot - 1013.0 + 1628.0 / x

  return
end
function f_09 ( x )

!*****************************************************************************80
!
!! f_09() evaluates function 9.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) f_09
  real ( kind = rk8 ) x

  f_09 = exp ( x ) - 2.0 - 0.01D+00 / x**2 + 0.000002D+00 / x**3

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
