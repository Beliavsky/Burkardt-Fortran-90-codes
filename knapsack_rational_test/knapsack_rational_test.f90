program main

!*****************************************************************************80
!
!! knapsack_rational_test() tests knapsack_rational().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 November 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'knapsack_rational_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test knapsack_rational().'

  call knapsack_rational_test ( )
  call knapsack_reorder_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'knapsack_rational_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' )  ''
  call timestamp ( )

  stop 0
end
subroutine knapsack_rational_test ( )

!*****************************************************************************80
!
!! knapsack_rational_test() tests knapsack_rational().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  integer i
  real ( kind = rk8 ) :: mass
  real ( kind = rk8 ) :: mass_limit = 26.0
  real ( kind = rk8 ), dimension ( n ) :: p = (/ &
    24.0, 13.0, 23.0, 15.0, 16.0 /)
  real ( kind = rk8 ) :: profit
  real ( kind = rk8 ), dimension ( n ) :: w = (/ &
    12.0,  7.0, 11.0,  8.0,  9.0 /)
  real ( kind = rk8 ), dimension ( n ) :: x

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'knapsack_rational_test():'
  write ( *, '(a)' ) '  knapsack_rational() solves the rational knapsack problem.'

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Object, Profit, Mass, "Profit Density"'
  write ( *, '(a)' ) ''
  do i = 1, n
    write ( *, '(i6,3f7.3)' ) i, p(i), w(i), p(i) / w(i)
  end do

  call knapsack_reorder ( n, p, w )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  After reordering by Profit Density:'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Object, Profit, Mass, "Profit Density"'
  write ( *, '(a)' ) ''
  do i = 1, n
    write ( *, '(i6,3f7.3)' ) i, p(i), w(i), p(i) / w(i)
  end do

  write ( *, '(a)' ) ''
  write ( *, '(a,f7.3)' ) '  Total mass restriction is ', mass_limit

  call knapsack_rational ( n, mass_limit, p, w, x, mass, profit )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Object, Density, Choice, Profit, Mass'
  write ( *, '(a)' ) ''
  do i = 1, n
    write ( *, '(i6,f7.3,f7.3,2f7.3)' ) i, p(i) / w(i), x(i), &
      x(i) * p(i), x(i) * w(i)
  end do

  write ( *, '(a)' ) ''
  write ( *, '(a,2f7.3)' ) '  Total:            ', profit, mass

  return
end
subroutine knapsack_reorder_test ( )

!*****************************************************************************80
!
!! knapsack_reorder_test() tests knapsack_reorder().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  integer i
  real ( kind = rk8 ), dimension ( n ) :: p = (/ &
    24.0, 13.0, 23.0, 15.0, 16.0 /)
  real ( kind = rk8 ), dimension ( n ) :: w = (/ &
    12.0,  7.0, 11.0,  8.0,  9.0 /)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'knapsack_reorder_test():'
  write ( *, '(a)' ) '  knapsack_reorder() reorders the knapsack data.'

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Object, Profit, Mass, "Profit Density"'
  write ( *, '(a)' ) ''
  do i = 1, n
    write ( *, '(2x,i6,2x,f7.3,2x,f7.3,2x,f7.3)' ) i, p(i), w(i), p(i)/w(i)
  end do

  call knapsack_reorder ( n, p, w )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  After reordering by Profit Density:'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Object, Profit, Mass, "Profit Density"'
  write ( *, '(a)' ) ''
  do i = 1, n
    write ( *, '(2x,i6,2x,f7.3,2x,f7.3,2x,f7.3)' ) i, p(i), w(i), p(i) / w(i)
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
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 ) time
  integer values(8)
  integer y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

