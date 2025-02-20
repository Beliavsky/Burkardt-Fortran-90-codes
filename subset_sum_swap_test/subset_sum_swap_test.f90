program main

!*****************************************************************************80
!
!! subset_sum_swap_test() tests subset_sum_swap().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'subset_sum_swap_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test subset_sum_swap().'

  call subset_sum_swap_tests ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'subset_sum_swap_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine subset_sum_swap_tests ( )

!*****************************************************************************80
!
!! subset_sum_swap_tests() tests subset_sum_swap().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer n
  integer t
  integer test
  integer, parameter :: test_num = 7
  integer, allocatable :: w(:)

  do test = 1, test_num

    if ( test == 1 ) then
      n = 8
      allocate ( w(1:n) )
      w = (/ 15, 22, 14, 26, 32, 9, 16, 8 /)
      t = 53
    else if ( test == 2 ) then
      n = 10
      allocate ( w(1:n) )
      w = (/ 267,  493,  869,  961, 1000, 1153, 1246, 1598, 1766, 1922 /)
      t = 5842
    else if ( test == 3 ) then
      n = 21
      allocate ( w(1:n) )
      w = (/  518533, 1037066, 2074132, 1648264, 796528, &
             1593056,  686112, 1372224,  244448, 488896, &
              977792, 1955584, 1411168,  322336, 644672, &
             1289344,   78688,  157376,  314752, 629504, &
             1259008 /)
      t = 2463098
    else if ( test == 4 ) then
      n = 10
      allocate ( w(1:n) )
      w = (/ 41, 34, 21, 20,  8,  7,  7,  4,  3,  3 /)
      t = 50
    else if ( test == 5 ) then
      n = 9
      allocate ( w(1:n) )
      w = (/ 81, 80, 43, 40, 30, 26, 12, 11, 9 /)
      t = 100
    else if ( test == 6 ) then
      n = 6
      allocate ( w(1:n) )
      w = (/ 1, 2, 4, 8, 16, 32 /)
      t = 22
    else if ( test == 7 ) then
      n = 10
      allocate ( w(1:n) )
      w = (/ 25, 27, 3, 12, 6, 15, 9, 30, 21, 19 /)
      t = 50
    end if

    call subset_sum_swap_try ( n, w, t )

    deallocate ( w )

  end do

  return
end
subroutine subset_sum_swap_try ( n, w, t )

!*****************************************************************************80
!
!! subset_sum_swap_try() tries the swap code for a given subset_sum problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer n

  integer i
  integer index(n)
  integer sum_achieved
  integer t
  integer w(n)

  write ( *, '(a)' ) ''
  write ( *, '(a,i8)' ) '  Target value: ', t
  
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Available weights:'
  write ( *, '(a)' ) ''
  do i = 1, n
    write ( *, '(a,i8)' ) '  ', w(i)
  end do

  call subset_sum_swap ( n, w, t, index, sum_achieved )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Selected weights:'
  write ( *, '(a)' ) ''

  do i = 1, n
    if ( index(i) == 1 ) then
      write ( *, '(a,i8)' ) '  ', w(i)
    end if
  end do

  write ( *, '(a)' ) ''
  write ( *, '(a,i8)' ) '  The target was      ', t
  write ( *, '(a,i8)' ) '  The achieved sum is ', sum_achieved
  write ( *, '(a,i8)' ) '  The defect is       ', t - sum_achieved

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

