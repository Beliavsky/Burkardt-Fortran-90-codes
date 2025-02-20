program main

!*****************************************************************************80
!
!! glomin_test() tests glomin().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 June 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) e
  real ( kind = rk ), external :: h_01
  real ( kind = rk ), external :: h_02
  real ( kind = rk ), external :: h_03
  real ( kind = rk ), external :: h_04
  real ( kind = rk ), external :: h_05
  real ( kind = rk ) m
  real ( kind = rk ) machep
  real ( kind = rk ) t

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'glomin_test():'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  glomin() seeks a global minimizer of a function F(X)'
  write ( *, '(a)' ) '  in an interval [A,B],'
  write ( *, '(a)' ) '  given some upper bound M for F".'

  machep = epsilon ( machep )
  e = sqrt ( machep )
  t = sqrt ( machep )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Tolerances'
  write ( *, '(a,g14.6)' ) '  e = ', e
  write ( *, '(a,g14.6)' ) '  t = ', t

  a = 7.0D+00
  b = 9.0D+00
  c = ( a + b ) / 2.0D+00
  m = 0.0D+00

  call glomin_example ( a, b, c, m, machep, e, t, h_01, &
    'h_01(x) = 2 - x' )

  a = 7.0D+00
  b = 9.0D+00
  c = ( a + b ) / 2.0D+00
  m = 100.0D+00

  call glomin_example ( a, b, c, m, machep, e, t, h_01, &
    'h_01(x) = 2 - x' )

  a = -1.0D+00
  b = +2.0D+00
  c = ( a + b ) / 2.0D+00
  m = 2.0D+00

  call glomin_example ( a, b, c, m, machep, e, t, h_02, &
    'h_02(x) = x * x' )

  a = -1.0D+00
  b = +2.0D+00
  c = ( a + b ) / 2.0D+00
  m = 2.1D+00

  call glomin_example ( a, b, c, m, machep, e, t, h_02, &
    'h_02(x) = x * x' )

  a = -0.5D+00
  b =  +2.0D+00
  c = ( a + b ) / 2.0D+00
  m = 14.0D+00

  call glomin_example ( a, b, c, m, machep, e, t, h_03, &
    'h_03(x) = x^3 + x^2' )

  a = -0.5D+00
  b =  +2.0D+00
  c = ( a + b ) / 2.0D+00
  m = 28.0D+00

  call glomin_example ( a, b, c, m, machep, e, t, h_03, &
    'h_03(x) = x^3 + x^2' )

  a = -10.0D+00
  b = +10.0D+00
  c = ( a + b ) / 2.0D+00
  m = 72.0D+00

  call glomin_example ( a, b, c, m, machep, e, t, h_04, &
    'h_04(x) = ( x + sin(x) ) * exp(-x*x)' )

  a = -10.0D+00
  b = +10.0D+00
  c = ( a + b ) / 2.0D+00
  m = 72.0D+00

  call glomin_example ( a, b, c, m, machep, e, t, h_05, &
    'h_05(x) = ( x - sin(x) ) * exp(-x*x)' )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'glomin_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine glomin_example ( a, b, c, m, machep, e, t, f, title )

!*****************************************************************************80
!
!! glomin_example() tests glomin() on one test function.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 May 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) A, B, the endpoints of the interval.
!
!    real ( kind = rk ) C, an initial guess for the global
!    minimizer.  If no good guess is known, C = A or B is acceptable.
!
!    real ( kind = rk ) M, the bound on the second derivative.
!
!    real ( kind = rk ) MACHEP, an estimate for the relative machine
!    precision.
!
!    real ( kind = rk ) E, a positive tolerance, a bound for the
!    absolute error in the evaluation of F(X) for any X in [A,B].
!
!    real ( kind = rk ) T, a positive absolute error tolerance.
!
!    external real ( kind = rk ) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose global minimum is being sought.
!
!    character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  integer calls
  real ( kind = rk ) e
  real ( kind = rk ), external :: f
  real ( kind = rk ) fa
  real ( kind = rk ) fb
  real ( kind = rk ) fx
  real ( kind = rk ) glomin
  real ( kind = rk ) m
  real ( kind = rk ) machep
  real ( kind = rk ) t
  character ( len = * ) title
  real ( kind = rk ) x

  fx = glomin ( a, b, c, m, machep, e, t, f, x, calls )
  fa = f ( a )
  fb = f ( b )

  write ( *, '(a)' ) ' '
  write ( *, '(2x,a)' ) trim ( title )
  write ( *, '(a,g14.6)' ) '  M = ', m
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      A                 X             B'
  write ( *, '(a)' ) '    F(A)              F(X)          F(B)'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,f14.8,2x,f14.8,2x,f14.8)' ) a,  x,  b
  write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) fa, fx, fb
  write ( *, '(a,i8)' ) '  Number of calls to F = ', calls

  return
end
function h_01 ( x )

!*****************************************************************************80
!
!! h_01() evaluates 2 - x.
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
!    real ( kind = rk ) H_01, the value of the function at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) h_01
  real ( kind = rk ) x

  h_01 = 2.0D+00 - x

  return
end
function h_02 ( x )

!*****************************************************************************80
!
!! h_02() evaluates x^2.
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
!    real ( kind = rk ) H_02, the value of the function at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) h_02
  real ( kind = rk ) x

  h_02 = x * x

  return
end
function h_03 ( x )

!*****************************************************************************80
!
!! h_03() evaluates x^3+x^2.
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
!    real ( kind = rk ) H_03, the value of the function at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) h_03
  real ( kind = rk ) x

  h_03 = x * x * ( x + 1.0D+00 )

  return
end
function h_04 ( x )

!*****************************************************************************80
!
!! h_04() evaluates ( x + sin ( x ) ) * exp ( - x * x ).
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
!    real ( kind = rk ) H_04, the value of the function at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) h_04
  real ( kind = rk ) x

  h_04 = ( x + sin ( x ) ) * exp ( - x * x )

  return
end
function h_05 ( x )

!*****************************************************************************80
!
!! h_05() evaluates ( x - sin ( x ) ) * exp ( - x * x ).
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
!    real ( kind = rk ) H_05, the value of the function at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) h_05
  real ( kind = rk ) x

  h_05 = ( x - sin ( x ) ) * exp ( - x * x )

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

