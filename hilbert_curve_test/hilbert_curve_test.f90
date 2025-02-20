program main

!*****************************************************************************80
!
!! hilbert_curve_test() tests hilbert_curve().
!
!  Modified:
!
!    02 January 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hilbert_curve_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test hilbert_curve().'

  call d2xy_test ( )
  call rot_test ( )
  call xy2d_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'HILBERT_CURVE_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine d2xy_test ( )

!*****************************************************************************80
!
!! D2XY_TEST tests D2XY.
!
!  Modified:
!
!    02 January 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer d
  integer m
  integer n
  integer x
  integer y

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'D2XY_TEST():'
  write ( *, '(a)' ) '  D2XY converts a Hilbert linear D coordinate to an (X,Y) 2D coordinate.'

  m = 3
  n = 2 ** m

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '    D    X    Y'
  write ( *, '(a)' ) ''
  do d = 0, n * n - 1
    call d2xy ( m, d, x, y )
    write ( *, '(2x,i3,2x,i3,2x,i3)' ) d, x, y
  end do

  return
end
subroutine rot_test ( )

!*****************************************************************************80
!
!! ROT_TEST tests ROT.
!
!  Modified:
!
!    02 January 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer m
  integer n
  integer rx
  integer ry
  integer x
  integer x0
  integer x1
  integer y
  integer y0
  integer y1

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ROT_TEST():'
  write ( *, '(a)' ) '  ROT rotates and flips a quadrant appropriately.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   X   Y  X0  Y0  X1  Y1'
  write ( *, '(a)' ) ''

  m = 3
  n = 2 ** m
  ry = 0

  do y = 0, 7
    do x = 0, 7
      rx = 0
      x0 = x
      y0 = y
      call rot ( n, x0, y0, rx, ry )
      rx = 1
      x1 = x
      y1 = y
      call rot ( n, x1, y1, rx, ry )
      write ( *, '(2x,i2,2x,i2,2x,i2,2x,i2,2x,i2,2x,i2)' ) x, y, x0, y0, x1, y1
    end do
  end do

  return
end
subroutine xy2d_test ( )

!*****************************************************************************80
!
!! XY2D_TEST tests XY2D.
!
!  Modified:
!
!    01 January 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer d
  integer m
  integer n
  integer x
  integer y

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'XY2D_TEST():'
  write ( *, '(a)' ) '  XY2D converts an (X,Y) 2D coordinate to a Hilbert linear D coordinate.'

  m = 3
  n = 2 ** m

  write ( *, '(a)' ) ''
  write ( *, '(a)', advance = 'no' ) '        '
  do x = 0, n - 1
    write ( *, '(i3)', advance = 'no' ) x
  end do
  write ( *, '(a)' ) ''

  write ( *, '(a)' ) ''
  do y = n - 1, 0, -1
    write ( *, '(2x,i3,a)', advance = 'no' ) y, ':  '
    do x = 0, n - 1
      call xy2d ( m, x, y, d )
      write ( *, '(i3)', advance = 'no' ) d
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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
