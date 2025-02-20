program main

!*****************************************************************************80
!
!! asa159_test() tests asa159().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa159_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test asa159().'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa159_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests RCONT2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 March 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 5
  integer, parameter :: n = 5

  integer a(m,n)
  integer, save, dimension ( n ) :: c = (/ 2, 2, 2, 2, 1 /)
  integer i
  integer ierror
  logical key
  integer, parameter :: ntest = 10
  integer, save, dimension ( m ) :: r = (/ 3, 2, 2, 1, 1 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  RCONT2 constructs a random matrix with'
  write ( *, '(a)' ) '  given row and column sums.'

  call i4vec_print ( m, r, '  The rowsum vector:' )
  call i4vec_print ( n, c, '  The columnsum vector:' )

  key = .false.

  do i = 1, ntest

    call rcont2 ( m, n, r, c, key, a, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  RCONT2 returned error flag IERROR = ', ierror
      return
    end if

    call i4mat_print ( m, n, a, '  The rowcolsum matrix:' )

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

