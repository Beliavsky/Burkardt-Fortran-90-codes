program main

!*****************************************************************************80
!
!! closest_pair_brute_test() tests closest_pair_brute().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 April 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer i
  integer, dimension ( 8 ) :: nvec = &
    (/ 100, 200, 400, 800, 1600, 3200, 6400, 12800 /)
  real ( kind = rk8 ) svec(8)

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'closest_pair_brute_test()'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test closest_pair_brute().'
!
!  Solve a sequence of problems of increasing size N.
!
  do i = 1, 8
    call closest_pair_brute_test01 ( nvec(i), svec(i) )
  end do
!
!  Print the final table.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Test  Points  Time'
  write ( *, '(a)' ) ''
  do i = 1, 8
    write ( *, '(2x,i4,2x,i6,2x,g14.6)' ) i, nvec(i), svec(i)
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'closest_pair_brute_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine closest_pair_brute_test01 ( n, s )

!*****************************************************************************80
!
!! closest_pair_brute_test01() tests the closest_pair_brute program.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n: the number of points.
!
!  Output:
!
!    real s: execution time in seconds.
!
  implicit none

  integer n

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer clock_count1
  integer clock_count2
  integer clock_max
  integer clock_rate
  real ( kind = rk8 ) d
  integer i
  integer j
  real ( kind = rk8 ) s
  real ( kind = rk8 ) xy(n,2)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'closest_pair_brute_test01():'
  write ( *, '(a)' ) '  closest_pair_brute() on a set of N points in 2D.'

  call random_number ( harvest = xy(1:n,1:2) )

  call system_clock ( clock_count1, clock_rate, clock_max )
  call closest_pair_brute ( n, xy, d, i, j )
  call system_clock ( clock_count2, clock_rate, clock_max )

  s = real ( clock_count2 - clock_count1, kind = rk8 ) / real ( clock_rate, kind = rk8 )

  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  n = ', n
  write ( *, '(a,g14.6)' ) '  time in seconds = ', s
  write ( *, '(a,i4,a,g14.6,g14.6)' ) '  xy(', i, ') = ', xy(i,1:2)
  write ( *, '(a,i4,a,g14.6,g14.6)' ) '  xy(', j, ') = ', xy(j,1:2)
  write ( *, '(a,g14.6)' ) '  distance = ', d

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
