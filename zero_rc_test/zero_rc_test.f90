program main

!*****************************************************************************80
!
!! zero_rc_test() tests zero_rc().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 December 2016
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
  real ( kind = rk ) t

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'zero_rc_test()()'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  zero_rc() seeks a root of a function F(X)'
  write ( *, '(a)' ) '  in an interval [A,B], using reverse communication.'

  t = epsilon ( t )

  a = 1.0D+00
  b = 2.0D+00
  call example_test ( a, b, t, f_01, 'f_01(x) = sin ( x ) - x / 2' )

  a = 0.0D+00
  b = 1.0D+00
  call example_test ( a, b, t, f_02, 'f_02(x) = 2 * x - exp ( - x )' )

  a = -1.0D+00
  b =  0.5D+00
  call example_test ( a, b, t, f_03, 'f_03(x) = x * exp ( - x )' )

  a =  0.0001D+00
  b =  20.0D+00
  call example_test ( a, b, t, f_04, 'f_04(x) = exp ( x ) - 1 / ( 100 * x * x )' )

  a = -5.0D+00
  b =  2.0D+00
  call example_test ( a, b, t, f_05, 'f_05(x) = (x+3) * (x-1) * (x-1)' )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'zero_rc_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine example_test ( a, b, t, f, title )

!*****************************************************************************80
!
!! example_test() tests zero() on one test function.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) A, B, the endpoints of the change of sign interval.
!
!    real ( kind = rk ) T, the error tolerance.
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
  real ( kind = rk ) arg
  real ( kind = rk ) b
  real ( kind = rk ), external :: f
  integer status
  real ( kind = rk ) t
  character ( len = *  ) title
  real ( kind = rk ) value
  
  write ( *, '(a)' ) ' '
  write ( *, '(2x,a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STATUS      X               F(X)'
  write ( *, '(a)' ) ' '

  status = 0

  do 

    call zero_rc ( a, b, t, arg, status, value )

    if ( status < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  zero_rc() returned an error flag!'
      exit
    end if

    value = f ( arg )

    write ( *, '(2x,i8,2x,g24.16,2x,g14.8)' ) status, arg, value

    if ( status == 0 ) then
      exit 
    end if

  end do

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

