program main

!*****************************************************************************80
!
!! doomsday_test() tests doomsday().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'doomsday_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test doomsday().'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DOOMSDAY_TEST():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests DOOMSDAY against a couple of test dates.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer d
  integer m
  integer w
  character ( len = 10 ) s1
  character ( len = 10 ) s2
  integer y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Try a couple selected dates.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YYYY  MM  DD  Weekday    Weekday'
  write ( *, '(a)' ) '                Tabulated  Computed'
  write ( *, '(a)' ) ' '

  y = 1989
  m = 7
  d = 13
  call doomsday_gregorian ( y, m, d, w )
  call weekday_to_name_common ( w, s1 )
  s2 = 'Thursday'
  write ( *, '(2x,i4,2x,i2,2x,i2,2x,a10,2x,a10)' ) y, m, d, s1, s2

  y = 2012
  m = 5
  d = 26
  call doomsday_gregorian ( y, m, d, w )
  call weekday_to_name_common ( w, s1 )
  s2 = 'Saturday'
  write ( *, '(2x,i4,2x,i2,2x,i2,2x,a10,2x,a10)' ) y, m, d, s1, s2

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests DOOMSDAY against a number of known values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer d
  integer m
  integer n_data
  integer w1
  integer w2
  character ( len = 10 ) s1
  character ( len = 10 ) s2
  integer y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  WEEKDAY_VALUES supplies a list of dates and weekdays.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YYYY  MM  DD  Weekday    Weekday'
  write ( *, '(a)' ) '                Tabulated  Computed'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call weekday_values ( n_data, y, m, d, w1 )

    if ( n_data <= 0 ) then
      exit
    end if
!
!  The transition from Julian to Gregorian calendars occurred in 1582
!  (for some people).  The data in "WEEKDAY_VALUES" before the transition
!  is stored in Julian format, which DOOMSDAY_GREGORIAN can't handle.
!  So let's just refuse to handle 1582 or earlier!
!
    if ( y <= 1582 ) then
      cycle
    end if

    call doomsday_gregorian ( y, m, d, w2 )

    call weekday_to_name_common ( w1, s1 )
    call weekday_to_name_common ( w2, s2 )

    write ( *, '(2x,i4,2x,i2,2x,i2,2x,a10,2x,a10)' ) y, m, d, s1, s2

  end do

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

