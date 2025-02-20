program main

!*****************************************************************************80
!
!! asa144_test() tests asa144().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa144_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test asa144().'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa144_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() tests rcont().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nrow = 5
  integer ( kind = 4 ), parameter :: ncol = 5

  integer ( kind = 4 ) ifault
  logical key
  integer ( kind = 4 ) matrix(nrow,ncol)
  integer ( kind = 4 ), dimension ( ncol ) :: ncolt = (/ &
    2, 2, 2, 2, 1 /)
  integer ( kind = 4 ), dimension ( nrow ) :: nrowt = (/ &
    3, 2, 2, 1, 1 /)
  integer ( kind = 4 ) nsubt(ncol)
  integer ( kind = 4 ) test
  integer ( kind = 4 ) :: test_num = 10

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  rcont() constructs a random matrix with'
  write ( *, '(a)' ) '  given row and column sums.'

  call i4vec_print ( nrow, nrowt, '  The rowsum vector:' )
  call i4vec_print ( ncol, ncolt, '  The columnsum vector: ' )

  key = .false.

  do test = 1, test_num

    call rcont ( nrow, ncol, nrowt, ncolt, nsubt, matrix, key, ifault )

    if ( ifault /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  RCONT returned IFAULT = ', ifault
      return
    end if

    call i4mat_print ( nrow, ncol, matrix, '  The rowcolsum matrix:' )

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

