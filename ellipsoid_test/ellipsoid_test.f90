program main

!*****************************************************************************80!
!
!! ellipsoid_test() tests ellipsoid().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 January 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ellipsoid_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test ellipsoid().'

  call ellipsoid_area_test ( )
  call ellipsoid_volume_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ellipsoid_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine ellipsoid_area_test ( )

!*****************************************************************************80
!
!! ellipsoid_area_test() tests ellipsoid_area().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 January 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) area
  real ( kind = rk ) ellipsoid_area

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ellipsoid_area_test():'
  write ( *, '(a)' ) '  ellipsoid_area() computes the surface area of an ellipsoid'
  write ( *, '(a)' ) '  satisfying the equation:'
  write ( *, '(a)' ) '    (x/a)^2+(y/b)^2+(z/c)^2=1'

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           A           B           C        Area'
  write ( *, '(a)' ) ''

  a = 1.0D+00
  b = 0.8D+00
  c = 0.625D+00
  area = ellipsoid_area ( a, b, c )
  write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) a, b, c, area

  a = 1.0D+00
  b = 1.0D+00
  c = 0.5D+00
  area = ellipsoid_area ( a, b, c )
  write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) a, b, c, area

  a = 1.0D+00
  b = 1.0D+00
  c = 1.0D+00
  area = ellipsoid_area ( a, b, c )
  write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) a, b, c, area

  a = 2.0D+00
  b = 1.0D+00
  c = 0.25D+00
  area = ellipsoid_area ( a, b, c )
  write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) a, b, c, area

  a = 2.0D+00
  b = 3.0D+00
  c = 4.0D+00
  area = ellipsoid_area ( a, b, c )
  write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) a, b, c, area

  return
end
subroutine ellipsoid_volume_test ( )

!*****************************************************************************80
!
!! ellipsoid_volume_test() tests ellipsoid_volume().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 January 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) ellipsoid_volume
  real ( kind = rk ) volume

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ellipsoid_volume_test():'
  write ( *, '(a)' ) '  ellipsoid_volume() computes the volume of an ellipsoid'
  write ( *, '(a)' ) '  satisfying the equation:'
  write ( *, '(a)' ) '    (x/a)^2+(y/b)^2+(z/c)^2=1'

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           A           B           C        Volume'
  write ( *, '(a)' ) ''

  a = 1.0D+00
  b = 0.8D+00
  c = 0.625D+00
  volume = ellipsoid_volume ( a, b, c )
  write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) a, b, c, volume

  a = 1.0D+00
  b = 1.0D+00
  c = 0.5D+00
  volume = ellipsoid_volume ( a, b, c )
  write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) a, b, c, volume

  a = 1.0D+00
  b = 1.0D+00
  c = 1.0D+00
  volume = ellipsoid_volume ( a, b, c )
  write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) a, b, c, volume

  a = 2.0D+00
  b = 1.0D+00
  c = 0.25D+00
  volume = ellipsoid_volume ( a, b, c )
  write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) a, b, c, volume

  a = 2.0D+00
  b = 3.0D+00
  c = 4.0D+00
  volume = ellipsoid_volume ( a, b, c )
  write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) a, b, c, volume

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

