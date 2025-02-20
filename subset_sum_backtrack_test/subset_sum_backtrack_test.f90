program main

!*****************************************************************************80
!
!! subset_sum_test() tests subset_sum().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 November 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'subset_sum_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test subset_sum().'

  call subset_sum_backtrack_tests ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'subset_sum_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine subset_sum_backtrack_tests ( )

!*****************************************************************************80
!
!! subset_sum_backtrack_tests() tests subset_sum_backtrack_test().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 November 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer n
  integer s
  integer, allocatable :: v(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'subset_sum_backtrack_tests():'
  write ( *, '(a)' ) '  subset_sum_backtrack_test() solves the subset sum problem'
  write ( *, '(a)' ) '  for specific values of S, N and V.'
  
  s = 9
  n = 5
  allocate ( v(1:n) )
  v = (/ 1, 2, 3, 5, 7 /)
  call subset_sum_backtrack_test ( s, n, v )
  deallocate ( v )
  
  s = 8
  n = 9
  allocate ( v(1:n) )
  v = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9 /)
  call subset_sum_backtrack_test ( s, n, v )
  deallocate ( v )
!
!  What happens with a repeated target?
!
  s = 8
  n = 9
  allocate ( v(1:n) )
  v = (/ 1, 2, 3, 3, 5, 6, 7, 8, 9 /)
  call subset_sum_backtrack_test ( s, n, v )
  deallocate ( v )
!
!  What happens with a target that needs all the values?
!
  s = 18
  n = 5
  allocate ( v(1:n) )
  v = (/ 1, 2, 3, 5, 7 /)
  call subset_sum_backtrack_test ( s, n, v )
  deallocate ( v )
!
!  A larger S.
!
  s = 5842
  n = 10
  allocate ( v(1:n) )
  v = (/ 267, 493, 869, 961, 1000, 1153, 1246, 1598, 1766, 1922 /)
  call subset_sum_backtrack_test ( s, n, v )
  deallocate ( v )

  return
end
subroutine subset_sum_backtrack_test ( s, n, v )

!*****************************************************************************80
!
!! subset_sum_backtrack_test() tests subset_sum_backtrack().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 November 2022
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer S, the desired sum.
!
!    Input, integer N, the number of values.
!
!    Input, integer V(N), the values.
!    These must be nonnegative, and sorted in ascending order.  
!    Duplicate values are allowed.
!
  implicit none

  integer n

  integer i
  integer k
  logical more
  logical plus
  integer s
  integer t
  integer u(n)
  integer v(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'subset_sum_backtrack_test():'
  write ( *, '(a)' ) '  subset_sum_backtrack() finds the "next" subset of the values V'
  write ( *, '(a)' ) '  which sum to the desired total S.'

  more = .false.
  u(1:n) = 0
  t = 0
  
  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  Desired sum S = ', s
  write ( *, '(a,i3)' ) '  Number of targets = ', n
  write ( *, '(a)', advance = 'no' ) '  Targets:'
  do i = 1, n
    write ( *, '(1x,i6)', advance = 'no' ) v(i)
  end do
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) ''

  k = 0
  
  do
    call subset_sum_backtrack ( s, n, v, more, u, t )
    if ( .not. more ) then
      exit
    end if
    k = k + 1
    write ( *, '(2x,i3,a,2x,i6,a)', advance = 'no' ) k, ':', s, ' = '
    plus = .false.
    do i = 1, n
      if ( u(i) /= 0 ) then
        if ( plus ) then
          write ( *, '(a)', advance = 'no' ) ' +'
        end if
        write ( *, '(i6)', advance = 'no' ) v(i)
        plus = .true.
      end if
    end do
    write ( *, '(a)' ) ''
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

