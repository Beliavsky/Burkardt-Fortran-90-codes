program main

!*****************************************************************************80
!
!! polygon_test() tests polygon().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 December 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'polygon_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test polygon().'

  call polygon_angles_test ( )
  call polygon_area_test ( )
  call polygon_area_lattice_test ( )
  call polygon_centroid_test ( )
  call polygon_contains_point_test ( )
  call polygon_convex_contains_point_test ( )
  call polygon_contains_point_3_test ( )
  call polygon_data_test ( )
  call polygon_diameter_test ( )
  call polygon_expand_test ( )
  call polygon_integral_test ( )
  call polygon_is_ccw_test ( )
  call polygon_is_convex_test ( )
  call polygon_perimeter_test ( )
  call polygon_perimeter_quad_test ( )
  call polygon_point_dist_test ( )
  call polygon_point_near_test ( )
  call polygon_sample_test ( )
  call polygon_triangulate_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'polygon_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
function f1 ( x, y )

!*****************************************************************************80
!
!! f1() evaluates f(x,y) = 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f1
  real ( kind = rk ) x
  real ( kind = rk ) y

  call r8_fake_use ( x )
  call r8_fake_use ( y )

  f1 = 1.0D+00

  return
end
function fx2 ( x, y )

!*****************************************************************************80
!
!! fx2() evaluates f(x,y) = x^2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fx2
  real ( kind = rk ) x
  real ( kind = rk ) y

  call r8_fake_use ( y )

  fx2 = x ** 2

  return
end
subroutine polygon_angles_test ( )

!*****************************************************************************80
!
!! polygon_angles_test() tests polygon_angles().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 6

  real ( kind = rk ) angle(n)
  integer i
  real ( kind = rk ) r8_degrees
  real ( kind = rk ), dimension (2,n) :: v = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    2.0D+00, 1.0D+00, &
    3.0D+00, 0.0D+00, &
    3.0D+00, 2.0D+00, &
    1.0D+00, 2.0D+00 /), (/ 2, n /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'polygon_angles_test():'
  write ( *, '(a)' ) '  polygon_angles() computes the angles of a polygon.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of polygonal vertices = ', n

  call r8mat_transpose_print ( 2, n, v, '  The polygon vertices:' )

  call polygon_angles ( n, v, angle )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Polygonal angles in degrees:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,2x,g14.6)' ) i, r8_degrees ( angle(i) )
  end do

  return
end
subroutine polygon_area_test ( )

!*****************************************************************************80
!
!! polygon_area_test() tests polygon_area().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 April 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: test_num = 2

  real ( kind = rk ) area1
  real ( kind = rk ) area2
  real ( kind = rk ) area_exact
  real ( kind = rk ), dimension ( test_num ) :: area_exact_test = (/ &
    2.0D+00, 6.0D+00 /)
  integer n
  integer, dimension ( test_num ) :: n_test = (/ 4, 8 /)
  integer test
  real ( kind = rk ), allocatable, dimension ( :, : ) :: v

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'polygon_area_test():'
  write ( *, '(a)' ) '  polygon_area()   computes the area of a polygon.'
  write ( *, '(a)' ) '  polygon_area_2() computes the area of a polygon.'

  do test = 1, test_num

    n = n_test(test)
    area_exact = area_exact_test(test)

    allocate ( v(1:2,1:n) )

    if ( test == 1 ) then

      v(1:2,1:n) = reshape ( (/ &
        1.0D+00, 0.0D+00, &
        2.0D+00, 1.0D+00, &
        1.0D+00, 2.0D+00, &
        0.0D+00, 1.0D+00 /), (/ 2, n /) )

    else if ( test == 2 ) then

      v(1:2,1:n) = reshape ( (/ &
        0.0D+00, 0.0D+00, &
        3.0D+00, 0.0D+00, &
        3.0D+00, 3.0D+00, &
        2.0D+00, 3.0D+00, &
        2.0D+00, 1.0D+00, &
        1.0D+00, 1.0D+00, &
        1.0D+00, 2.0D+00, &
        0.0D+00, 2.0D+00 /), (/ 2, n /) )

    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of polygonal vertices = ', n

    call r8mat_transpose_print ( 2, n, v, '  The polygon vertices:' )

    call polygon_area ( n, v, area1 )
    call polygon_area_2 ( n, v, area2 )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Exact area is        ', area_exact
    write ( *, '(a,g14.6)' ) '  polygon_area():      ', area1
    write ( *, '(a,g14.6)' ) '  polygon_area_2():    ', area2
 
    deallocate ( v )

  end do

  return
end
subroutine polygon_area_lattice_test ( )

!*****************************************************************************80
!
!! polygon_area_lattice_test() tests polygon_area_lattice().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  integer b
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'polygon_area_lattice_test():'
  write ( *, '(a)' ) '  polygon_area_lattice() returns the area'
  write ( *, '(a)' ) '  of a polygon, measured in lattice points.'

  i = 5
  b = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Number of interior lattice points = ', i
  write ( *, '(a,i4)' ) '  Number of boundary lattice points = ', b

  call polygon_area_lattice ( i, b, area )

  write ( *, '(a,g14.6)' ) '  Area of polygon is ', area

  return
end
subroutine polygon_centroid_test ( )

!*****************************************************************************80
!
!! polygon_centroid_test() tests polygon_centroid().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk ) centroid1(2)
  real ( kind = rk ) centroid2(2)
  real ( kind = rk ), dimension (2,n) :: v = reshape ( (/ &
    1.0D+00, 0.0D+00, &
    2.0D+00, 1.0D+00, &
    1.0D+00, 2.0D+00, &
    0.0D+00, 1.0D+00 /), (/ 2, n /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'polygon_centroid_test():'
  write ( *, '(a)' ) '  polygon_centroid()   computes the centroid of a polygon.'
  write ( *, '(a)' ) '  polygon_centroid_2() computes the centroid of a polygon.'

  call r8mat_transpose_print ( 2, n, v, '  The polygon vertices:' )

  call polygon_centroid ( n, v, centroid1 )
  call r8vec_print ( 2, centroid1, '  polygon_centroid():' )

  call polygon_centroid ( n, v, centroid2 )
  call r8vec_print ( 2, centroid2, '  polygon_centroid_2():' )

  return
end
subroutine polygon_contains_point_test ( )

!*****************************************************************************80
!
!! polygon_contains_point_test() tests polygon_contains_point().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 April 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5
  integer, parameter :: test_num = 4

  logical inside
  real ( kind = rk ) p(2)
  real ( kind = rk ), dimension (2,test_num) :: p_test = reshape ( (/ &
    1.0D+00,  1.0D+00, &
    3.0D+00,  4.0D+00, &
    0.0D+00,  2.0D+00, &
    0.5D+00, -0.25D+00 /), (/2, test_num /) )
  integer test
  real ( kind = rk ), dimension(2,n) :: v = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    2.0D+00, 1.0D+00, &
    1.0D+00, 2.0D+00, &
    0.0D+00, 2.0D+00 /), (/ 2, n /) )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'polygon_contains_point_test():'
  write ( *, '(a)' ) '  polygon_contains_point()   determines if a point is in a polygon.'

  call r8mat_transpose_print ( 2, n, v, '  The polygon vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          P          Inside?'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
 
    p(1:2) = p_test(1:2,test)
 
    call polygon_contains_point ( n, v, p, inside )

    write ( *, '(2x,2g14.6,4x,l1)' ) p(1:2), inside

  end do
 
  return
end
subroutine polygon_convex_contains_point_test ( )

!*****************************************************************************80
!
!! polygon_convex_contains_point_test() tests polygon_convex_contains_point().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 January 2025
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5
  integer, parameter :: test_num = 4

  logical inside
  real ( kind = rk ) p(2)
  real ( kind = rk ), dimension (2,test_num) :: p_test = reshape ( (/ &
    1.0D+00,  1.0D+00, &
    3.0D+00,  4.0D+00, &
    0.0D+00,  2.0D+00, &
    0.5D+00, -0.25D+00 /), (/2, test_num /) )
  integer test
  real ( kind = rk ), dimension(2,n) :: v = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    2.0D+00, 1.0D+00, &
    1.0D+00, 2.0D+00, &
    0.0D+00, 2.0D+00 /), (/ 2, n /) )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'polygon_convex_contains_point_test():'
  write ( *, '(a)' ) '  polygon_convex_contains_point() determines if'
  write ( *, '(a)' ) '  a point is in a convex polygon.'

  call r8mat_transpose_print ( 2, n, v, '  The polygon vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          P          Inside?'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
 
    p(1:2) = p_test(1:2,test)
 
    call polygon_convex_contains_point ( n, v, p, inside )

    write ( *, '(2x,2g14.6,4x,l1)' ) p(1:2), inside

  end do
 
  return
end
subroutine polygon_contains_point_3_test ( )

!*****************************************************************************80
!
!! polygon_contains_point_3_test() tests polygon_contains_point_3().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 April 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5
  integer, parameter :: test_num = 4

  logical inside
  real ( kind = rk ) p(2)
  real ( kind = rk ), dimension (2,test_num) :: p_test = reshape ( (/ &
    1.0D+00,  1.0D+00, &
    3.0D+00,  4.0D+00, &
    0.0D+00,  2.0D+00, &
    0.5D+00, -0.25D+00 /), (/2, test_num /) )
  integer test
  real ( kind = rk ), dimension(2,n) :: v = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    2.0D+00, 1.0D+00, &
    1.0D+00, 2.0D+00, &
    0.0D+00, 2.0D+00 /), (/ 2, n /) )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'polygon_contains_point_3_test():'
  write ( *, '(a)' ) '  polygon_contains_point_3() determines if a point is in a polygon.'

  call r8mat_transpose_print ( 2, n, v, '  The polygon vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          P          Inside?'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
 
    p(1:2) = p_test(1:2,test)
 
    call polygon_contains_point_3 ( n, v, p, inside )

    write ( *, '(2x,2g14.6,4x,l1)' ) p(1:2), inside

  end do
 
  return
end
subroutine polygon_data_test ( )

!*****************************************************************************80
!
!! polygon_data_test() tests polygon_data_*();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 April 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  integer n
  real ( kind = rk ) radin
  real ( kind = rk ) radout
  real ( kind = rk ) side

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'polygon_data_test():'
  write ( *, '(a)' ) '  polygon_data_inrad() uses the inradius for polygon data.'
  write ( *, '(a)' ) '  polygon_data_outrad() uses the outradius for polygon data.'
  write ( *, '(a)' ) '  polygon_data_side() uses the side for polygon data.'

  do n = 3, 5

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of polygonal sides = ', n

    side = 1.0D+00
    write ( *, '(a)' ) ' '

    call polygon_data_side ( n, side, area, radin, radout )
    write ( *, '(a)' ) ''
    write ( *, '(a,g14.6)' )   '    AREA =   ', area
    write ( *, '(a,g14.6)' )   '    RADIN =  ', radin
    write ( *, '(a,g14.6)' )   '    RADOUT = ', radout
    write ( *, '(a,g14.6,a)' ) '    SIDE =   ', side, '  (given)'

    call polygon_data_inrad ( n, radin, area, radout, side )
    write ( *, '(a)' ) ''
    write ( *, '(a,g14.6)' )   '    AREA =   ', area
    write ( *, '(a,g14.6,a)' ) '    RADIN =  ', radin, '  (given)'
    write ( *, '(a,g14.6)' )   '    RADOUT = ', radout
    write ( *, '(a,g14.6)' )   '    SIDE =   ', side

    call polygon_data_outrad ( n, radout, area, radin, side )
    write ( *, '(a)' ) ''
    write ( *, '(a,g14.6)' )   '    AREA =   ', area
    write ( *, '(a,g14.6)' )   '    RADIN =  ', radin
    write ( *, '(a,g14.6,a)' ) '    RADOUT = ', radout, '  (given)'
    write ( *, '(a,g14.6)' )   '    SIDE =   ', side

  end do

  return
end
subroutine polygon_diameter_test ( )

!*****************************************************************************80
!
!! polygon_diameter_test() tests polygon_diameter();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk ) diameter
  real ( kind = rk ) :: diameter_exact = 2.0D+00
  real ( kind = rk ), dimension ( 2, n ) :: v = reshape ( (/ &
    1.0D+00, 0.0D+00, &
    2.0D+00, 1.0D+00, &
    1.0D+00, 2.0D+00, &
    0.0D+00, 1.0D+00 /), (/ 2, n /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'polygon_diameter_test():'
  write ( *, '(a)' ) '  polygon_diameter() computes the diameter of a polygon.'

  call r8mat_transpose_print ( 2, n, v, '  The polygon vertices:' )

  call polygon_diameter ( n, v, diameter )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Diameter ( computed ) ', diameter
  write ( *, '(a,g14.6)' ) '  Diameter ( exact )    ', diameter_exact
 
  return
end
subroutine polygon_expand_test ( )

!*****************************************************************************80
!
!! polygon_expand_test() tests polygon_expand();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk ) h
  real ( kind = rk ), dimension(2,n) :: v = reshape ( (/ &
    1.0D+00, 1.0D+00, &
    5.0D+00, 1.0D+00, &
    2.0D+00, 4.0D+00, &
    1.0D+00, 3.0D+00 /), (/ 2, n /) )
  real ( kind = rk ), dimension(2,n) :: w

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'polygon_expand_test():'
  write ( *, '(a)' ) '  polygon_expand() "expands" a polygon by an amount H.'

  h = 0.5D+00

  call r8mat_transpose_print ( 2, n, v, '  The polygon vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The expansion amount H = ', h

  call polygon_expand ( n, v, h, w )

  call r8mat_transpose_print ( 2, n, w, '  The expanded polygon:' )

  return
end
subroutine polygon_integral_test ( )

!*****************************************************************************80
!
!! polygon_integral_test() tests polygon_integral_*().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 April 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3
 
  real ( kind = rk ) result
  real ( kind = rk ), dimension(2,n) :: v = reshape ( (/ &
    1.0D+00, 1.0D+00, &
    4.0D+00, 3.0D+00, &
    2.0D+00, 5.0D+00 /), (/ 2, n /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'polygon_integral_test():'
  write ( *, '(a)' ) '  polygon_integral_1()  integrates 1 over a polygon.'
  write ( *, '(a)' ) '  polygon_integral_x()  integrates x over a polygon.'
  write ( *, '(a)' ) '  polygon_integral_y()  integrates y over a polygon.'
  write ( *, '(a)' ) '  polygon_integral_xx() integrates xx over a polygon.'
  write ( *, '(a)' ) '  polygon_integral_xy() integrates xy over a polygon.'
  write ( *, '(a)' ) '  polygon_integral_yy() integrates yy over a polygon.'

  call r8mat_transpose_print ( 2, n, v, '  Vertices of polygon V:' )

  write ( *, '(a)' ) ''
  call polygon_integral_1 ( n, v, result )
  write ( *, '(a,4x,g14.6)' ) '  Integral of 1  over V', result

  call polygon_integral_x ( n, v, result )
  write ( *, '(a,4x,g14.6)' ) '  Integral of X  over V', result

  call polygon_integral_y ( n, v, result )
  write ( *, '(a,4x,g14.6)' ) '  Integral of Y  over V', result

  call polygon_integral_xx ( n, v, result )
  write ( *, '(a,4x,g14.6)' ) '  Integral of XX over V', result

  call polygon_integral_xy ( n, v, result )
  write ( *, '(a,4x,g14.6)' ) '  Integral of XY over V', result

  call polygon_integral_yy ( n, v, result )
  write ( *, '(a,4x,g14.6)' ) '  Integral of YY over V', result

  return
end
subroutine polygon_is_ccw_test ( )

!*****************************************************************************80
!
!! polygon_is_ccw_test() tests polygon_is_ccw().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 December 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk ) a(2)
  real ( kind = rk ) b(2)
  real ( kind = rk ) c(2)
  real ( kind = rk ) d(2)
  integer test
  real ( kind = rk ) v(n,2)
  logical value

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'polygon_is_ccw_test():'
  write ( *, '(a)' ) '  polygon_is_cc2() determines if the vertices of a'
  write ( *, '(a)' ) '  polygon are listed in counter-clockwise order.'

  a = (/ 0.0, 0.0 /)
  b = (/ 1.0, 0.0 /)
  c = (/ 1.0, 1.0 /)
  d = (/ 0.0, 1.0 /)
!
!  Try all six orientations of the the four vertices.
!  Only one (a,b,c,d) should be acceptable.
!
  do test = 1, 6
 
    if ( test == 1 ) then
      v(1,1:2) = a
      v(2,1:2) = b
      v(3,1:2) = c
      v(4,1:2) = d
    else if ( test == 2 ) then
      v(1,1:2) = a
      v(2,1:2) = b
      v(3,1:2) = d
      v(4,1:2) = c
    else if ( test == 3 ) then
      v(1,1:2) = a
      v(2,1:2) = c
      v(3,1:2) = b
      v(4,1:2) = d
    else if ( test == 4 ) then
      v(1,1:2) = a
      v(2,1:2) = c
      v(3,1:2) = d
      v(4,1:2) = b
    else if ( test == 5 ) then
      v(1,1:2) = a
      v(2,1:2) = d
      v(3,1:2) = b
      v(4,1:2) = c
    else if ( test == 6 ) then
      v(1,1:2) = a
      v(2,1:2) = d
      v(3,1:2) = c
      v(4,1:2) = b
    end if

    call r8mat_print ( n, 2, v, '  Polygon vertices:' )
 
    call polygon_is_ccw ( n, v, value )

    if ( value ) then
      write ( *, '(a)' ) '  The polygon vertices are counter clockwise.'
    else
      write ( *, '(a)' ) '  The polygon vertices are NOT counter clockwise.'
    end if

  end do

  return
end
subroutine polygon_is_convex_test ( )

!*****************************************************************************80
!
!! polygon_is_convex_test() tests polygon_is_convex().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 10
  integer, parameter :: test_num = 11
 
  real ( kind = rk ) angle
  integer i
  character ( len = 80 ) message(-1:2)
  integer n
  integer polygon_is_convex
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  integer result
  integer test
  character ( len = 255 ) title
  real ( kind = rk ) v(2,n_max)

  message(-1) = 'The polygon is not convex.'
  message( 0) = 'The polygon is degenerate and convex.'
  message( 1) = 'The polygon is convex and counterclockwise.'
  message( 2) = 'The polygon is convex and clockwise.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'polygon_is_convex_test():'
  write ( *, '(a)' ) '  polygon_is_convex() determines if a polygon is convex.'

  do test = 1, test_num

    if ( test == 1 ) then
      n = 1
      v(1:2,1:n) = reshape ( (/ &
        0.0D+00, 0.0D+00 /), (/ 2, n /) )
      title = '  A point:'
    else if ( test == 2 ) then
      n = 2
      v(1:2,1:n) = reshape ( (/ &
        0.0D+00, 0.0D+00, &
        1.0D+00, 2.0D+00 /), (/ 2, n /) )
      title = '  A line:'
    else if ( test == 3 ) then
      n = 3
      v(1:2,1:n) = reshape ( (/ &
        0.0D+00, 0.0D+00, &
        2.0D+00, 0.0D+00, &
        1.0D+00, 0.0D+00 /), (/ 2, n /) )
      title = '  A triangle:'
    else if ( test == 4 ) then
      n = 3
      v(1:2,1:n) = reshape ( (/ &
        0.0D+00, 0.0D+00, &
        1.0D+00, 0.0D+00, &
        0.0D+00, 2.0D+00 /), (/ 2, n /) )
      title = '  A CCW triangle:'
    else if ( test == 5 ) then
      n = 3
      v(1:2,1:n) = reshape ( (/ &
        0.0D+00, 0.0D+00, &
        0.0D+00, 2.0D+00, &
        1.0D+00, 0.0D+00 /), (/ 2, n /) )
      title = '  A CW triangle:'
    else if ( test == 6 ) then
      n = 4
      v(1:2,1:n) = reshape ( (/ &
        1.0D+00, 0.0D+00, &
        2.0D+00, 0.0D+00, &
        3.0D+00, 1.0D+00, &
        0.0D+00, 1.0D+00 /), (/ 2, n /) )
      title = '  Polygon with large angle:'
    else if ( test == 7 ) then
      n = 5
      v(1:2,1:n) = reshape ( (/ &
        0.0D+00, 0.0D+00, &
        0.5D+00, 0.5D+00, &
        1.0D+00, 0.0D+00, &
        1.0D+00, 1.0D+00, &
        0.0D+00, 1.0D+00 /), (/ 2, n /) )
      title = '  Polygon with huge angle:'
    else if ( test == 8 ) then
      n = 5
      do i = 1, n
        angle = real ( i - 1, kind = rk ) * 4.0D+00 * r8_pi &
          / real ( n, kind = rk )
        v(1,i) = cos ( angle )
        v(2,i) = sin ( angle )
      end do
      title = '  A five-pointed star:'
    else if ( test == 9 ) then
      n = 6
      do i = 1, n
        angle = real ( i - 1, kind = rk ) * 2.0D+00 * r8_pi &
          / real ( n, kind = rk )
        v(1,i) = cos ( angle )
        v(2,i) = sin ( angle )
      end do
      title = '  A hexagon:'
    else if ( test == 10 ) then
      n = 6
      v(1:2,1:n) = reshape ( (/ &
        0.0D+00, 0.0D+00, &
        2.0D+00, 0.0D+00, &
        1.0D+00, 1.0D+00, &
        0.0D+00, 0.0D+00, &
        2.0D+00, 0.0D+00, &
        1.0D+00, 1.0D+00 /), (/ 2, n /) )
      title = '  A triangle twice:'
    else if ( test == 11 ) then
      n = 8
      v(1:2,1:n) = reshape ( (/ &
        1.0D+00, 0.0D+00, &
        3.0D+00, 0.0D+00, &
        3.0D+00, 3.0D+00, &
        0.0D+00, 3.0D+00, &
        0.0D+00, 1.0D+00, &
        2.0D+00, 1.0D+00, &
        2.0D+00, 2.0D+00, &
        1.0D+00, 2.0D+00 /), (/ 2, n /) )
      title = '  Square knot:'
    end if

    call r8mat_transpose_print ( 2, n, v, title )
    result = polygon_is_convex ( n, v )
    write ( *, '(2x,a)' ) message(result)

  end do

  return
end
subroutine polygon_perimeter_test ( )

!*****************************************************************************80
!
!! polygon_perimeter_test() tests polygon_perimeter().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 4
  integer, parameter :: n2 = 3
 
  real ( kind = rk ) perimeter
  real ( kind = rk ), dimension(2,n1) :: v1 = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, &
    0.0D+00, 1.0D+00 /), (/ 2, n1 /) )
  real ( kind = rk ), dimension(2,n2) :: v2 = reshape ( (/ &
    1.0D+00, 1.0D+00, &
    4.0D+00, 3.0D+00, &
    2.0D+00, 5.0D+00 /), (/ 2, n2 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'polygon_perimeter_test():'
  write ( *, '(a)' ) '  polygon_perimeter() computes the perimeter of a polygon.'

  call r8mat_transpose_print ( 2, n1, v1, '  Vertices of polygon V1:' )

  call polygon_perimeter ( n1, v1, perimeter )
  write ( *, '(a)' ) ''
  write ( *, '(a,4x,g14.6)' ) '  Perimeter of V1 = ', perimeter

  call r8mat_transpose_print ( 2, n2, v2, '  Vertices of polygon V2:' )

  call polygon_perimeter ( n2, v2, perimeter )
  write ( *, '(a)' ) ''
  write ( *, '(a,4x,g14.6)' ) '  Perimeter of V2 = ', perimeter

  return
end
subroutine polygon_perimeter_quad_test ( )

!*****************************************************************************80
!
!! polygon_perimeter_quad_test() tests polygon_perimeter_quad().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 4
  integer, parameter :: n2 = 3
 
  real ( kind = rk ), external :: f1
  real ( kind = rk ), external :: fx2
  real ( kind = rk ) hmax
  integer i
  real ( kind = rk ), dimension(2,n1) :: v1 = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, &
    0.0D+00, 1.0D+00 /), (/ 2, n1 /) )
  real ( kind = rk ), dimension(2,n2) :: v2 = reshape ( (/ &
    1.0D+00, 1.0D+00, &
    4.0D+00, 3.0D+00, &
    2.0D+00, 5.0D+00 /), (/ 2, n2 /) )
  real ( kind = rk ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'polygon_perimeter_quad_test():'
  write ( *, '(a)' ) '  polygon_perimeter_quad() estimates the integral of'
  write ( *, '(a)' ) '  a function over the perimeter of a polygon using'
  write ( *, '(a)' ) '  the composite midpoint rule over each side.'

  call r8mat_transpose_print ( 2, n1, v1, '  Vertices of polygon V1:' )

  hmax = 0.5D+00
  call polygon_perimeter_quad ( n1, v1, hmax, f1, value )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6,a,g14.6)' ) '  Using HMAX = ', hmax, ' estimated integral of 1 over perimeter = ', value

  write ( *, '(a)' ) ''
  hmax = 2.0D+00
  do i = 1, 3
    hmax = hmax / 2.0D+00
    call polygon_perimeter_quad ( n1, v1, hmax, fx2, value )
    write ( *, '(a,g14.6,a,g14.6)' ) '  Using HMAX = ', hmax, ' estimated integral of x^2 over perimeter = ', value
  end do

  call r8mat_transpose_print ( 2, n2, v2, '  Vertices of polygon V2:' )

  hmax = 0.5D+00
  call polygon_perimeter_quad ( n2, v2, hmax, f1, value )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6,a,g14.6)' ) '  Using HMAX = ', hmax, ' estimated integral of 1 over perimeter = ', value

  write ( *, '(a)' ) ''
  hmax = 2.0D+00
  do i = 1, 3
    hmax = hmax / 2.0D+00
    call polygon_perimeter_quad ( n2, v2, hmax, fx2, value )
    write ( *, '(a,g14.6,a,g14.6)' ) '  Using HMAX = ', hmax, ' estimated integral of x^2 over perimeter = ', value
  end do

  return
end
subroutine polygon_point_dist_test ( )

!*****************************************************************************80
!
!! polygon_point_dist_test() tests polygon_point_dist().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3
 
  real ( kind = rk ) dist
  real ( kind = rk ) p(2)
  real ( kind = rk ), dimension ( 2, 3 ) :: p_test = reshape ( (/ &
    4.0D+00,  5.0D+00, &
    2.0D+00,  3.0D+00, &
   -2.0D+00, -1.0D+00 /), (/ 2, 3 /) )
  real ( kind = rk ), dimension(2,n) :: v = reshape ( (/ &
    1.0D+00, 1.0D+00, &
    4.0D+00, 3.0D+00, &
    2.0D+00, 5.0D+00 /), (/ 2, n /) )
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POLYGON_POINT_DIST_test():'
  write ( *, '(a)' ) '  POLYGON_POINT_DIST() computes polygon-point distance.'

  call r8mat_transpose_print ( 2, n, v, '  Vertices of polygon:' )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '       X             Y             DIST'
  write ( *, '(a)' ) ''

  do test = 1, 3

    p(1:2) = p_test(1:2,test)
    call polygon_point_dist ( n, v, p, dist )
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) p(1:2), dist

  end do

  return
end
subroutine polygon_point_near_test ( )

!*****************************************************************************80
!
!! polygon_point_near_test() tests polygon_point_near().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3
 
  real ( kind = rk ) dist
  real ( kind = rk ) p(2)
  real ( kind = rk ), dimension ( 2, 3 ) :: p_test = reshape ( (/ &
    4.0D+00,  5.0D+00, &
    2.0D+00,  3.0D+00, &
   -2.0D+00, -1.0D+00 /), (/ 2, 3 /) )
  real ( kind = rk ) pn(2)
  real ( kind = rk ), dimension(2,n) :: v = reshape ( (/ &
    1.0D+00, 1.0D+00, &
    4.0D+00, 3.0D+00, &
    2.0D+00, 5.0D+00 /), (/ 2, n /) )
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POLYGON_POINT_NEAR_test():'
  write ( *, '(a)' ) '  POLYGON_POINT_NEAR() computes nearest point on polygon.'

  call r8mat_transpose_print ( 2, n, v, '  Vertices of polygon:' )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '       X             Y             XN             YN'
  write ( *, '(a)' ) ''

  do test = 1, 3

    p(1:2) = p_test(1:2,test)
    call polygon_point_near ( n, v, p, pn, dist )
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) p(1:2), pn(1:2)

  end do

  return
end
subroutine polygon_sample_test ( )

!*****************************************************************************80
!
!! polygon_sample_test() tests polygon_sample().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 20
  integer, parameter :: nv = 6

  real ( kind = rk ), dimension(2,nv) :: v = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    2.0D+00, 0.0D+00, &
    2.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, &
    1.0D+00, 2.0D+00, &
    0.0D+00, 1.0D+00 /), (/ 2, nv /) )
  real ( kind = rk ) x(2,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POLYGON_SAMPLE_test():'
  write ( *, '(a)' ) '  POLYGON_SAMPLE() samples a polygon.'

  call polygon_sample ( nv, v, n, x )

  call r8mat_transpose_print ( 2, n, x, '  Sample points:' )

  return
end
subroutine polygon_triangulate_test ( )

!*****************************************************************************80
!
!! polygon_triangulate_test() tests polygon_triangulate().
!
!  Discussion:
!
!    There are N-3 triangles in the triangulation.
!
!    For the first N-2 triangles, the first edge listed is always an
!    internal diagonal.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 October 2015
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  integer j
  integer triangles(3,n-2)
  real ( kind = rk ), dimension ( n ) :: x = (/ &
    8.0D+00, 7.0D+00, 6.0D+00, 5.0D+00, 4.0D+00, &
    3.0D+00, 2.0D+00, 1.0D+00, 0.0D+00, 4.0D+00 /)
  real ( kind = rk ), dimension ( n ) :: y = (/ &
    0.0D+00, 10.0D+00,  0.0D+00, 10.0D+00,  0.0D+00, &
   10.0D+00,  0.0D+00, 10.0D+00,  0.0D+00, -2.0D+00 /)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'polygon_triangulate_test():'
  write ( *, '(a)' ) '  polygon_triangulate() triangulates a polygon.'
  write ( *, '(a)' ) '  Here, we triangulate the comb_10 polygon.'

  call polygon_triangulate ( n, x, y, triangles )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Triangles:'
  write ( *, '(a)' ) ''
  do j = 1, n - 2
    write ( *, '(2x,i2,4x,i2,2x,i2,2x,i2)' ) j, triangles(1:3,j)
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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
