program main

!*****************************************************************************80
!
!! levenshtein_matrix_test() tests levenshtein_matrix().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 September 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, allocatable :: d(:,:)
  integer i
  integer j
  integer m
  integer n
  character ( len = 5 ) :: s1 = 'water'
  character ( len = 6 ) :: s2 = 'kitten'
  character ( len = 8 ) :: s3 = 'saturday'
  character ( len = 10 ) :: s4 = 'pheromones'
  character ( len = 4 ) :: t1 = 'wine'
  character ( len = 7 ) :: t2 = 'sitting'
  character ( len = 6 ) :: t3 = 'sunday'
  character ( len = 12 ) :: t4 = 'photographer'

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'levenshtein_matrix_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  levenshtein_matrix() computes the Levenshtein matrix'
  write ( *, '(a)' ) '  associated with the computation of the Levenshtein'
  write ( *, '(a)' ) '  distance between two strings.'

  m = len ( s1 )
  n = len ( t1 )
  allocate ( d(0:m,0:n) )
  call levenshtein_matrix ( m, s1, n, t1, d )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  S = "' // s1 // '"'
  write ( *, '(a)' ) '  T = "' // t1 // '"'
  do i = 0, m
    do j = 0, n
      write ( *, '(1x,i2)', advance='no' ) d(i,j)
    end do
    write ( *, '(a)' ) ''
  end do
  deallocate ( d )

  m = len ( s2 )
  n = len ( t2 )
  allocate ( d(0:m,0:n) )
  call levenshtein_matrix ( m, s2, n, t2, d )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  S = "' // s2 // '"'
  write ( *, '(a)' ) '  T = "' // t2 // '"'
  do i = 0, m
    do j = 0, n
      write ( *, '(1x,i2)', advance='no' ) d(i,j)
    end do
    write ( *, '(a)' ) ''
  end do
  deallocate ( d )

  m = len ( s3 )
  n = len ( t3 )
  allocate ( d(0:m,0:n) )
  call levenshtein_matrix ( m, s3, n, t3, d )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  S = "' // s3 // '"'
  write ( *, '(a)' ) '  T = "' // t3 // '"'
  do i = 0, m
    do j = 0, n
      write ( *, '(1x,i2)', advance='no' ) d(i,j)
    end do
    write ( *, '(a)' ) ''
  end do
  deallocate ( d )

  m = len ( s4 )
  n = len ( t4 )
  allocate ( d(0:m,0:n) )
  call levenshtein_matrix ( m, s4, n, t4, d )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  S = "' // s4 // '"'
  write ( *, '(a)' ) '  T = "' // t4 // '"'
  do i = 0, m
    do j = 0, n
      write ( *, '(1x,i2)', advance='no' ) d(i,j)
    end do
    write ( *, '(a)' ) ''
  end do
  deallocate ( d )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'levenshtein_matrix_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
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
!    18 May 2013
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
