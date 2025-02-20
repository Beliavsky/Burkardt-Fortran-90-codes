program main

!*****************************************************************************80
!
!! middle_square_test() tests middle_square().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 September 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'middle_square_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test middle_square()'

  call middle_square_next_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'middle_square_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine middle_square_next_test ( )

!*****************************************************************************80
!
!! middle_square_next_test() tests middle_square_next().
!
!  Discussion:
!
!    This function simply demonstrates the results of 10 successive calls
!    to middle_square_next(), for a range of values of d.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 September 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: i20 = selected_int_kind ( 20 )
  integer ( kind = i20 ) d
  integer ( kind = i20 ) i
  integer ( kind = i20 ) s

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'middle_square_next_test():'
  write ( *, '(a)' ) '  middle_square_next ( s, d ) applies the middle square algorithm'
  write ( *, '(a)' ) '  using a 2*d digit seed.'

  do d = 1, 5

    write ( *, '(a)' )
    write ( *, '(a,i2)' ) '  Using d = ', d
    write ( *, '(a)' )

    s = mod ( 2147483647, 10**(2*d) )

    do i = 0, 10
      write ( *, '(2x,i2,2x,i12)' ) i, s
      call middle_square_next ( s, d, s )
    end do

  end do

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

