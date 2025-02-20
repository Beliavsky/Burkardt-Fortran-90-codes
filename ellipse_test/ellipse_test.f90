program main

!*****************************************************************************80
!
!! ellipse_test() tests ellipse().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 October 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ellipse_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test ellipse().'

  call ellipse_area1_test ( )
  call ellipse_area2_test ( )
  call ellipse_area3_test ( )
  call ellipse_aspect_ratio_test ( )
  call ellipse_eccentricity_test ( )
  call ellipse_flattening_test ( )
  call ellipse_point_dist_2d_test ( )
  call ellipse_point_near_2d_test ( )
  call ellipse_points_2d_test ( )
  call ellipse_points_arc_2d_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ellipse_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine ellipse_area1_test ( )

!*****************************************************************************80
!
!! ellipse_area1_test() tests ellipse_area1().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 April 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(2,2)
  real ( kind = rk ) area
  real ( kind = rk ) ellipse_area1
  integer i
  real ( kind = rk ) r

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ellipse_area1_test():'
  write ( *, '(a)' ) '  ellipse_area1() computes the area of an ellipse.'

  r = 10.0D+00

  a = reshape ( (/ 5.0D+00, 1.0D+00, 1.0D+00, 2.0D+00 /), (/ 2, 2 /) )

  area = ellipse_area1 ( a, r )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  R = ', r
  write ( *, '(a)' ) '  Matrix A in ellipse definition x*A*x=r^2'
  do i = 1, 2
    write ( *, '(2x,g14.6,2x,g14.6)' ) a(i,1:2)
  end do
  write ( *, '(a,g14.6)' ) '  Area = ', area

  return
end
subroutine ellipse_area2_test ( )

!*****************************************************************************80
!
!! ellipse_area2_test() tests ellipse_area2().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 April 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) area
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) ellipse_area2
 
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ellipse_area2_test():'
  write ( *, '(a)' ) '  ellipse_area2() computes the area of an ellipse.'

  a = 5.0D+00
  b = 2.0D+00
  c = 2.0D+00
  d = 10.0D+00

  area = ellipse_area2 ( a, b, c, d )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6,a,g14.6,a,g14.6,a,g14.6)' ) &
    '  Ellipse: ', a, ' * x^2 + ', b, ' * xy + ', c, ' * y^2 = ', d
  write ( *, '(a,g14.6)' ) '  Area = ', area

  return
end
subroutine ellipse_area3_test ( )

!*****************************************************************************80
!
!! ellipse_area3_test() tests ellipse_area3().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 April 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ) ellipse_area3
  real ( kind = rk ) r1
  real ( kind = rk ) r2
 
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ellipse_area3_test():'
  write ( *, '(a)' ) '  ellipse_area3() computes the area of an ellipse.'

  r1 = 10.0D+00
  r2 = 10.0D+00 / 3.0D+00

  area = ellipse_area3 ( r1, r2 )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6,a,g14.6,a)' ) &
    '  Ellipse: x^2/', r1, '^2 + y^2 / ', r2, '^2 = 1'
  write ( *, '(a,g14.6)' ) '  Area = ', area

  return
end
subroutine ellipse_aspect_ratio_test ( )

!*****************************************************************************80
!
!! ellipse_aspect_ratio_test() tests ellipse_aspect_ratio().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 October 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) e
  real ( kind = rk ) ellipse_aspect_ratio
  integer i
  integer n
 
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ellipse_aspect_ratio_test():'
  write ( *, '(a)' ) '  ellipse_aspect_ratio() computes the aspect ratio of an ellipse.'

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '      A      B      Ratio'
  write ( *, '(a)' ) ''
  a = 1.0D+00
  n = 10
  do i = 0, n
    b = real ( i, kind = rk ) / real ( n, kind = rk )
    e = ellipse_aspect_ratio ( a, b )
    write ( *, '(  f5.1  f5.1  f10.6)' ) a, b, e
  end do

  return
end
subroutine ellipse_eccentricity_test ( )

!*****************************************************************************80
!
!! ellipse_eccentricity_test() tests ellipse_eccentricity().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 April 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) e
  real ( kind = rk ) ellipse_eccentricity
  integer i
  integer n
 
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ellipse_eccentricity_test():'
  write ( *, '(a)' ) '  ellipse_eccentricity() computes the eccentricity of an ellipse.'

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '      A      B      Ecc'
  write ( *, '(a)' ) ''
  a = 1.0D+00
  n = 10
  do i = 0, n
    b = real ( i, kind = rk ) / real ( n, kind = rk )
    e = ellipse_eccentricity ( a, b )
    write ( *, '(  f5.1  f5.1  f10.6)' ) a, b, e
  end do

  return
end
subroutine ellipse_flattening_test ( )

!*****************************************************************************80
!
!! ellipse_flattening_test() tests ellipse_flattening().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 October 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) e
  real ( kind = rk ) ellipse_flattening
  integer i
  integer n
 
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ellipse_flattening_test():'
  write ( *, '(a)' ) '  ellipse_flattening() computes the flattening of an ellipse.'

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '      A      B      Flat'
  write ( *, '(a)' ) ''
  a = 1.0D+00
  n = 10
  do i = 0, n
    b = real ( i, kind = rk ) / real ( n, kind = rk )
    e = ellipse_flattening ( a, b )
    write ( *, '(  f5.1  f5.1  f10.6)' ) a, b, e
  end do

  return
end
subroutine ellipse_point_dist_2d_test ( )

!*****************************************************************************80
!
!! ellipse_point_dist_2d_test() tests ellipse_point_dist_2d().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 April 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) dist
  integer i
  integer, parameter :: n = 10
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) :: r1 = 3.0D+00
  real ( kind = rk ) :: r2 = 2.0D+00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ellipse_point_dist_2d_test():'
  write ( *, '(a)' ) '  ellipse_point_dist_2d() is given a point P, and'
  write ( *, '(a)' ) '  finds the distance to an ellipse in 2D.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  The ellipse is (X/R1)^2 + (Y/R2)^2 = 1'
  write ( *, '(a)' ) ''
  write ( *, '(a,f14.6)' ) '  R1 = ', r1
  write ( *, '(a,f14.6)' ) '  R2 = ', r2
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           P            DIST'
  write ( *, '(a)' ) ''

  do i = -3, n + 3

    p(1) = ( real ( n - i, kind = rk ) * 0.0D+00   &
           + real (     i, kind = rk ) * 4.0D+00 ) &
           / real ( n,     kind = rk )

    p(2) = ( real ( n - i, kind = rk ) * 3.0D+00   &
           + real (     i, kind = rk ) * 0.0D+00 ) &
           / real ( n,     kind = rk )

    call ellipse_point_dist_2d ( r1, r2, p, dist )

    write ( *, '(2x,2f8.4,2x,f8.4)' ) p(1:dim_num), dist

  end do

  return
end
subroutine ellipse_point_near_2d_test ( )

!*****************************************************************************80
!
!! ellipse_point_near_2d_test() tests ellipse_point_near_2d().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 April 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  integer i
  integer, parameter :: n = 10
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) :: r1 = 3.0D+00
  real ( kind = rk ) :: r2 = 2.0D+00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ellipse_point_near_2d_test():'
  write ( *, '(a)' ) '  ellipse_point_near_2d() is given a point P, and'
  write ( *, '(a)' ) '  finds the nearest point PN on an ellipse in 2D.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  The ellipse is (X/R1)^2 + (Y/R2)^2 = 1'
  write ( *, '(a)' ) ''
  write ( *, '(a,f14.6)' ) '  R1 = ', r1
  write ( *, '(a,f14.6)' ) '  R2 = ', r2
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           P                PN'
  write ( *, '(a)' ) ''

  do i = -3, n + 3

    p(1) = ( real ( n - i, kind = rk ) * 0.0D+00   &
           + real (     i, kind = rk ) * 4.0D+00 ) &
           / real ( n,     kind = rk )

    p(2) = ( real ( n - i, kind = rk ) * 3.0D+00   &
           + real (     i, kind = rk ) * 0.0D+00 ) &
           / real ( n,     kind = rk )

    call ellipse_point_near_2d ( r1, r2, p, pn )

    write ( *, '(2x,2f8.4,2x,2f8.4)' ) p(1:dim_num), pn(1:dim_num)

  end do

  return
end
subroutine ellipse_points_2d_test ( )

!*****************************************************************************80
!
!! ellipse_points_2d_test() tests ellipse_points_2d();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 April 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer, parameter :: n_max = 24

  real ( kind = rk ) area
  real ( kind = rk ) ellipse_area3
  integer j
  integer n
  real ( kind = rk ), dimension ( dim_num ) :: pc = (/ 5.0D+00, -2.0D+00 /)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) psi
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) v(dim_num,n_max)

  r1 = 3.0D+00
  r2 = 1.0D+00
  psi = r8_pi / 6.0D+00
  n = 16

  area = ellipse_area3 ( r1, r2 )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ellipse_points_2d_test():'
  write ( *, '(a)' ) '  ellipse_points_2d() returns points on an ellipse;'

  write ( *, '(a,g14.6,a,g14.6)' ) '  Ellipse center at ', pc(1), ', ', pc(2)

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6,a,g14.6)' ) '  radii R1 = ', r1, ' R2 = ', r2
  write ( *, '(a,g14.6)' ) '  and angle PSI = ', psi
  write ( *, '(a,g14.6)' ) '  and area = ', area

  call ellipse_points_2d ( pc, r1, r2, psi, n, v )

  do j = 1, n
    write ( *, '(2x,g14.6,2x,g14.6)' ) v(1:2,j)
  end do
 
  return
end
subroutine ellipse_points_arc_2d_test ( )

!*****************************************************************************80
!
!! ellipse_points_arc_2d_test() tests ellipse_points_arc_2d().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 April 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 13
  integer, parameter :: dim_num = 2

  integer j
  real ( kind = rk ) p(dim_num,n)
  real ( kind = rk ), dimension(dim_num) :: pc = (/ 5.0D+00, -2.0D+00 /)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) psi
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2

  r1 = 3.0D+00
  r2 = 1.0D+00
  psi = r8_pi / 6.0D+00
  theta1 = r8_pi / 2.0D+00
  theta2 = 2.0D+00 * r8_pi

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ellipse_points_arc_2d_test():'
  write ( *, '(a)' ) '  ellipse_points_arc_2d() returns points on an'
  write ( *, '(a)' ) '  elliptical arc.'
  write ( *, '(a)' ) ''
  write ( *, '(a,2g14.6)' ) '  The ellipse has center ', pc(1), pc(2)
  write ( *, '(a,g14.6,a,g14.6)' ) '  radii R1 = ', r1, ' R2 = ', r2
  write ( *, '(a,g14.6)' ) '  and angle PSI = ', psi
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  The arc extends from THETA1 = ', theta1
  write ( *, '(a,g14.6)') '  to THETA2 = ', theta2

  call ellipse_points_arc_2d ( pc, r1, r2, psi, theta1, theta2, n, p )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Sample points:'
  write ( *, '(a)' ) ''
  do j = 1, n
    write ( *, '(2x,g14.6,2x,g14.6)' ) p(1:2,j)
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

