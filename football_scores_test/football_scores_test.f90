program main

!*****************************************************************************80
!
!! football_scores_test() tests football_scores().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 September 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: i16 = selected_int_kind ( 16 )

  integer, parameter :: n_max = 50

  integer n
  integer ( kind = i16 ) s(0:n_max)

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'football_scores_test()'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  football_scores() computes the number of ways'
  write ( *, '(a)' ) '  of achieving a particular score in football.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We assume scoring options are:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  +1 for a one point safety (returned conversion'
  write ( *, '(a)' ) '     after the other team scores a touchdown.'
  write ( *, '(a)' ) '  +2 for a safety;'
  write ( *, '(a)' ) '  +3 for a field goal;'
  write ( *, '(a)' ) '  +6 for a touchdown with no followup;'
  write ( *, '(a)' ) '  +7 for a touchdown with a point bonus;'
  write ( *, '(a)' ) '  +8 for a touchdown with two point conversion;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Score                  Ways'
  write ( *, '(a)' ) ' '
!
!  You get a score of N points by
!    N-1 + a 1 point safety,
!    N-2 + a safety, or
!    N-3 + a field goal, or
!    N-6 + a touchdown with no followup
!    N-7 + a touchdown with a point bonus, or
!    N-8 + a touchdown with two point conversion.
!
  call football_scores ( n_max, s )

  do n = 0, n_max
    write ( *, '(2x,i6,2x,i25)' ) n, s(n)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'football_scores_test()'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
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

  integer, parameter :: rk = kind ( 1.0D+00 )

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

