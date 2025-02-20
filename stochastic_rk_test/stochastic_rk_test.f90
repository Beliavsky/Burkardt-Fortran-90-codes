program main

!*****************************************************************************80
!
!! Stochastic_rk_test() tests stochastic_rk().
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'stochastic_rk_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  test stochastic_rk().'
 
  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'stochastic_rk_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() tests rk1_ti_step().
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk8 ), external :: fi
  real ( kind = rk8 ), external :: gi
  real ( kind = rk8 ) h
  integer i
  real ( kind = rk8 ) q
  real ( kind = rk8 ) t
  real ( kind = rk8 ), parameter :: t0 = 0.0D+00
  real ( kind = rk8 ), parameter :: tn = 1.0D+00
  real ( kind = rk8 ) x(0:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  rk1_ti_step() uses a first order RK method'
  write ( *, '(a)' ) '  for a problem whose right hand side does not'
  write ( *, '(a)' ) '  depend explicitly on time.'

  h = ( tn - t0 ) / real ( n, kind = rk8 )
  q = 1.0D+00

  i = 0
  t = t0
  x(i) = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I           T             X'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i8,2x,f14.6,2x,g14.6)' ) i, t, x(i)

  do i = 1, n

    t = ( real ( n - i, kind = rk8 ) * t0   &
        + real (     i, kind = rk8 ) * tn ) &
        / real ( n,     kind = rk8 )

    call rk1_ti_step ( x(i-1), t, h, q, fi, gi, x(i) )

    write ( *, '(2x,i8,2x,f14.6,2x,g14.6)' ) i, t, x(i)

  end do

  return
end
function fi ( x )

!*****************************************************************************80
!
!! fi() is a time invariant deterministic right hand side.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk8 ) X, the argument.
!
!    Output, real ( kind = rk8 ) FI, the value.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fi
  real ( kind = rk8 ) x

  fi = 1.0D+00

  return
end
function gi ( x )

!*****************************************************************************80
!
!! gi() is a time invariant stochastic right hand side.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk8 ) X, the argument.
!
!    Output, real ( kind = rk8 ) GI, the value.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) gi
  real ( kind = rk8 ) x

  gi = 1.0D+00

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
!    This code is distributed under the GNU LGPL license.
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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

