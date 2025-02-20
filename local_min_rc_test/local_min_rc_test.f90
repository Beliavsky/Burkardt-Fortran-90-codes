program main

!*****************************************************************************80
!
!! local_min_rc_test() tests local_min_rc().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 May 2021
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
  write ( *, '(a)' ) 'local_min_rc_test():'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  local_min_rc() is a reverse communication code '
  write ( *, '(a)' ) '  which seeks a local minimizer of a function F(X)'
  write ( *, '(a)' ) '  in an interval [A,B].'

  a = 0.0D+00
  b = 3.141592653589793D+00
  call example_test ( a, b, g_01, 'g_01(x) = ( x - 2 ) * ( x - 2 ) + 1' )

  a = 0.0D+00
  b = 1.0D+00
  call example_test ( a, b, g_02, 'g_02(x) = x * x + exp ( - x )' )

  a = -2.0D+00
  b =  2.0D+00
  call example_test ( a, b, g_03, 'g_03(x) = x^4 + 2x^2 + x + 3' )

  a =  0.0001D+00
  b =  1.0D+00
  call example_test ( a, b, g_04, 'g_04(x) = exp ( x ) + 1 / ( 100 x )' )

  a =  0.0002D+00
  b = 2.0D+00
  call example_test ( a, b, g_05, 'g_05(x) = exp ( x ) - 2x + 1/(100x) - 1/(1000000x^2)' )

  a = 1.8D+00
  b = 1.9D+00
  call example_test ( a, b, g_06, 'g_06(x) = - x * sin ( 10 pi x ) - 1' )

  a = 0.0D+00
  b = 2.0D+00
  call example_test ( a, b, g_07, 'g_07(x) = 2x^4 - 4x^2 + x + 20' )

  a = -2.0D+00
  b = 0.0D+00
  call example_test ( a, b, g_07, 'g_07(x) = 2x^4 - 4x^2 + x + 20' )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'local_min_rc_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine example_test ( a, b, f, title )

!*****************************************************************************80
!
!! example_test() tests local_min_rc() on one example.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 May 2021
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
  real ( kind = rk ) a2
  real ( kind = rk ) arg
  real ( kind = rk ) b
  real ( kind = rk ) b2
  real ( kind = rk ), external :: f
  integer ( kind = 4 ) status
  integer ( kind = 4 ) step
  real ( kind = rk ) t
  character ( len = * ) title
  real ( kind = rk ) value

  t = 10.0D+00 * sqrt ( epsilon ( t ) )

  write ( *, '(a)' ) ' '
  write ( *, '(2x,a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Step      X                          F(X)'
  write ( *, '(a)' ) ' '
  step = 0

  arg = a
  value = f ( arg )
  write ( *, '(2x,i4,2x,g24.16,2x,g24.16)' ) step, arg, value

  arg = b
  value = f ( arg )
  write ( *, '(2x,i4,2x,g24.16,2x,g24.16)' ) step, arg, value

  a2 = a
  b2 = b
  status = 0

  do

    call local_min_rc ( a2, b2, arg, status, value )
 
    if ( status < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'example_test(): Fatal error!'
      write ( *, '(a)' ) '  local_min_rc() returned negative status.'
      exit
    end if

    value = f ( arg )

    step = step + 1
    write ( *, '(2x,i4,2x,g24.16,2x,g24.16)' ) step, arg, value

    if ( status == 0 ) then
      exit
    end if 

  end do

  return
end
function g_01 ( x )

!*****************************************************************************80
!
!! G_01 evaluates (x-2)^2 + 1.
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
!  Parameters:
!
!    Input, real ( kind = rk ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = rk ) G_01, the value of the function at X.
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
!! G_02 evaluates x^2 + exp ( - x ).
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
!  Parameters:
!
!    Input, real ( kind = rk ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = rk ) G_02, the value of the function at X.
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
!! G_03 evaluates x^4+2x^2+x+3.
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
!  Parameters:
!
!    Input, real ( kind = rk ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = rk ) G_03, the value of the function at X.
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
!! G_04 evaluates exp(x)+1/(100X)
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
!  Parameters:
!
!    Input, real ( kind = rk ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = rk ) G_04, the value of the function at X.
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
!! G_05 evaluates exp(x) - 2x + 1/(100x) - 1/(1000000x^2)
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
!  Parameters:
!
!    Input, real ( kind = rk ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = rk ) G_05, the value of the function at X.
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
!! G_06 evaluates - x * sin ( 10 pi x ) - 1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 November 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = rk ) G_06, the value of the function at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) g_06
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) x

  g_06 = - x * sin ( 10.0 * r8_pi * x ) - 1.0D+00

  return
end
function g_07 ( x )

!*****************************************************************************80
!
!! G_07 evaluates 2x^4 - 4x^2 + x + 20
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 September 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = rk ) G_07, the value of the function at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) g_07
  real ( kind = rk ) x

  g_07 = 2.0D+00 * x**4 - 4.0D+00 * x**2 + x + 20.0D+00

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
!  Parameters:
!
!    None
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

