program main

!*****************************************************************************80
!
!! knapsack_values_test() tests knapsack_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 November 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer i
  integer k
  integer n
  integer n_data
  real ( kind = rk8 ) r
  integer, allocatable :: s(:)
  integer, allocatable :: v(:)
  integer, allocatable :: w(:)

  call timestamp ( )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'knapsack_values_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test knapsack_values().'

  n_data = 0

  do
!
!  First call returns value of n.
!
    n = 0
    call knapsack_values ( n_data, n, v, w, s, k )

    if ( n == 0 ) then
      exit
    end if

    allocate ( v(1:n) )
    allocate ( w(1:n) )
    allocate ( s(1:n) )
!
!  Second call returns data, and increments n_data.
!
    call knapsack_values ( n_data, n, v, w, s, k )

    write ( *, '(a)' ) ''
    write ( *, '(a,i4)' ) '  Problem #', n_data
    write ( *, '(a,i4)' ) '  Number of items is ', n
    write ( *, '(a,i6)' ) '  Knapsack weight limit is ', k
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) '   Item 0/1  Value  Weight  Value/Weight'
    write ( *, '(a)' ) ''
    do i = 1, n
      r = real ( v(i), kind = rk8 ) / real ( w(i), kind = rk8 )
      write ( *, '(2x,i5,2x,i2,2x,i8,2x,i8,2x,f7.2)' ) &
        i, s(i), v(i), w(i), r
    end do
    write ( *, '(a)' ) ''
    write ( *, '(a,2x,i2,2x,i8,2x,i8,2x,f7.2)' ) &
      '  Taken', &
      sum ( s(1:n) ), &
      dot_product ( s, v ), &
      dot_product ( s, w ), &
      real ( dot_product ( s, v ), kind = rk8 ) / &
      real ( dot_product ( s, w ), kind = rk8 )

    deallocate ( s )
    deallocate ( w )
    deallocate ( v )

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'knapsack_values_test():'
  write ( *, '(a)' ) '  Normal end of execution.'

  call timestamp ( )

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
