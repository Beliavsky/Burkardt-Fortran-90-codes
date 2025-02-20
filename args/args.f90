program main

!*****************************************************************************80
!
!! args() tests the (semi-standard) GETARGS utility.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Usage:
!
!    args arg1 arg2 arg3 ...
!
  implicit none

  character ( len = 255 ) arg
  integer i
  integer iargc
  integer numarg

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'args():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate the use of command line argument'
  write ( *, '(a)' ) '  routines in a FORTRAN program.'
  write ( *, '(a)' ) '  These include getarg() and iargc().'

  numarg = iargc ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) &
    '  ARGS was called with IARGC() = ', numarg, ' arguments.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CALL GETARG(I,ARG) returns the arguments:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I     ARG '
  write ( *, '(a)' ) ' '

  do i = 0, numarg
    call getarg ( i, arg )
    write ( *, '(2x,i3,2x,a20)' ) i, arg
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'args():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
!    26 August 2021
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
