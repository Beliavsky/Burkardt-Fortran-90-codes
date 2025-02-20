program main

!*****************************************************************************80
!
!! triangle_test() tests triangle().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 January 2025
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'triangle_test():'
  write ( *, '(a)' ) '  Fortran90 version:'
  write ( *, '(a)' ) '  Test triangle().'

  call triangle_angles_test ( )
  call triangle_area_test ( )
  call triangle_centroid_test ( )
  call triangle_circumcircle_test ( )
  call triangle_contains_point_test ( )
  call triangle_diameter_test ( )
  call triangle_edge_length_test ( )
  call triangle_incircle_test ( )
  call triangle_orientation_test ( )
  call triangle_orthocenter_test ( )
  call triangle_point_dist_test ( )
  call triangle_point_near_test ( )
  call triangle_quality_test ( )
  call triangle_reference_sample_test ( )
  call triangle_sample_test ( )
  call triangle_sample_reflect_test ( )
  call triangle_xsi_to_xy_test ( )
  call triangle_xy_to_xsi_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'triangle_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine triangle_angles_test ( )

!*****************************************************************************80
!
!! triangle_angles_test() tests triangle_angles().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) angle(3)
  integer i
  real ( kind = rk8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk8 ), dimension ( 2, 3 ) :: t = reshape ( (/ &
    0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00 /), (/ 2, 3 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'triangle_angles_test():'
  write ( *, '(a)' ) '  triangle_angles() computes the angles in a triangle.'

  call r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  call triangle_angles ( t, angle )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Radians      Degrees'
  write ( *, '(a)' ) ' '
  do i = 1, 3
    write ( *, '(2x,g14.6,2x,g14.6)' ) angle(i), ( angle(i) ) * 180.0D+00 / r8_pi
  end do

  return
end
subroutine triangle_area_test ( )

!*****************************************************************************80
!
!! triangle_area_test() tests triangle_area().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) area
  real ( kind = rk8 ), dimension ( 2, 3 ) :: t = reshape ( (/ &
    0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00 /), (/ 2, 3 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE_AREA_TEST():'
  write ( *, '(a)' ) '  TRIANGLE_AREA() computes the area of a triangle.'

  call r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  call triangle_area ( t, area )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Triangle area is ', area

  return
end
subroutine triangle_centroid_test ( )

!*****************************************************************************80
!
!! TRIANGLE_CENTROID_TEST() tests TRIANGLE_CENTROID();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: test_num = 4

  real ( kind = rk8 ) centroid(2)
  real ( kind = rk8 ) t(2,3)
  real ( kind = rk8 ), dimension(2,3,test_num) :: t_test = reshape ( (/ &
         0.0D+00,  0.0D+00, &
         1.0D+00,  0.0D+00, &
         0.0D+00,  1.0D+00, &
         0.0D+00,  0.0D+00, &
         1.0D+00,  0.0D+00, &
         0.5D+00,  0.86602539D+00, &
         0.0D+00,  0.0D+00, &
         1.0D+00,  0.0D+00, &
         0.5D+00, 10.0D+00, &
         0.0D+00,  0.0D+00, &
         1.0D+00,  0.0D+00, &
        10.0D+00,  2.0D+00 /), (/ 2, 3, test_num /) )
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE_CENTROID_TEST():'
  write ( *, '(a)' ) '  TRIANGLE_CENTROID() computes the centroid of a triangle.'

  do test = 1, test_num

    t(1:2,1:3) = t_test(1:2,1:3,test)

    call r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

    call triangle_centroid ( t, centroid )

    call r8vec_print ( 2, centroid, '  Centroid:' )
 
  end do

  return
end
subroutine triangle_circumcircle_test ( )

!*****************************************************************************80
!
!! TRIANGLE_CIRCUMCIRCLE_TEST() tests TRIANGLE_CIRCUMCIRCLE().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: test_num = 4

  real ( kind = rk8 ) pc(2)
  real ( kind = rk8 ) r
  real ( kind = rk8 ) t(2,3)
  real ( kind = rk8 ), dimension(2,3,test_num) :: t_test = reshape ( (/ &
         0.0D+00,  0.0D+00, &
         1.0D+00,  0.0D+00, &
         0.0D+00,  1.0D+00, &
         0.0D+00,  0.0D+00, &
         1.0D+00,  0.0D+00, &
         0.5D+00,  0.86602539D+00, &
         0.0D+00,  0.0D+00, &
         1.0D+00,  0.0D+00, &
         0.5D+00, 10.0D+00, &
         0.0D+00,  0.0D+00, &
         1.0D+00,  0.0D+00, &
        10.0D+00,  2.0D+00 /), (/ 2, 3, test_num /) )
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE_CIRCUMCIRCLE_TEST():'
  write ( *, '(a)' ) '  TRIANGLE_CIRCUMCIRCLE() computes the circumcircle of a triangle.'

  do test = 1, test_num

    t(1:2,1:3) = t_test(1:2,1:3,test)

    call r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

    call triangle_circumcircle ( t, r, pc )

    call r8vec_print ( 2, pc, '  Circumcenter' )

    write ( *, '(a)' ) ''
    write ( *, '(a,g14.6)' ) '  Circumradius: ', r

  end do

  return
end
subroutine triangle_contains_point_test ( )

!*****************************************************************************80
!
!! TRIANGLE_CONTAINS_POINT_TEST() tests TRIANGLE_CONTAINS_POINT().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: test_num = 7

  logical inside
  integer j
  real ( kind = rk8 ) p(2)
  real ( kind = rk8 ), dimension ( 2, test_num ) :: p_test = reshape ( (/ &
     0.25D+00,   0.25D+00, &
     0.75D+00,   0.25D+00, &
     1.00D+00,   1.00D+00, &
    11.00D+00,   0.50D+00, &
     0.00D+00,   1.00D+00, &
     0.50D+00, -10.00D+00, &
     0.60D+00,   0.60D+00 /), (/ 2, test_num /) )
  real ( kind = rk8 ), dimension ( 2, 3 ) :: t = reshape ( (/ &
    0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00 /), (/ 2, 3 /) )
  real ( kind = rk8 ), dimension ( 2, 3 ) :: t2
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE_CONTAINS_POINT_TEST():'
  write ( *, '(a)' ) '  TRIANGLE_CONTAINS_POINT() reports if a point '
  write ( *, '(a)' ) '  is inside a triangle.'

  call r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X       Y     Inside'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:2) = p_test(1:2,test)
 
    call triangle_contains_point ( t, p, inside )

    write ( *, '(2x,2f8.3,5x,l1)' ) p(1:2), inside

  end do
!
!  Make a copy of the triangle with vertices in reverse order.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeat the test, but reverse the triangle vertex'
  write ( *, '(a)' ) '  ordering.'
 
  do j = 1, 3
    t2(1:2,j) = t(1:2,4-j)
  end do

  call r8mat_transpose_print ( 2, 3, t2, &
    '  Triangle vertices (reversed):' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X       Y     Inside'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:2) = p_test(1:2,test)
 
    call triangle_contains_point ( t2, p, inside )

    write ( *, '(2x,2f8.3,5x,l1)' ) p(1:2), inside

  end do
 
  return
end
subroutine triangle_diameter_test ( )

!*****************************************************************************80
!
!! TRIANGLE_DIAMETER_TEST() tests TRIANGLE_DIAMETER().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: test_num = 3

  real ( kind = rk8 ) diameter
  real ( kind = rk8 ) t(2,3)
  real ( kind = rk8 ), dimension(2,3,test_num) :: t_test = reshape ( (/ &
     4.0D+00, 2.0D+00, &
     1.0D+00, 5.0D+00, &
    -2.0D+00, 2.0D+00, &
     4.0D+00, 2.0D+00, &
     5.0D+00, 4.0D+00, &
     6.0D+00, 6.0D+00, &
     4.0D+00, 2.0D+00, &
     1.0D+00, 5.0D+00, &
     4.0D+00, 2.0D+00 /), (/ 2, 3, test_num /) )
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE_DIAMETER_TEST():'
  write ( *, '(a)' ) '  TRIANGLE_DIAMETER() computes the diameter of '
  write ( *, '(a)' ) '  the SMALLEST circle around the triangle.'

  do test = 1, test_num

    t(1:2,1:3) = t_test(1:2,1:3,test)

    call r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

    call triangle_diameter ( t, diameter )

    write ( *, '(a)' ) ''
    write ( *, '(a,g14.6)' ) '  Diameter =       ', diameter

  end do

  return
end
subroutine triangle_edge_length_test ( )

!*****************************************************************************80
!
!! TRIANGLE_EDGE_LENGTH_TEST() tests TRIANGLE_EDGE_LENGTH().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: test_num = 3

  real ( kind = rk8 ) edge_length(3)
  real ( kind = rk8 ) t(2,3)
  real ( kind = rk8 ), dimension(2,3,test_num) :: t_test = reshape ( (/ &
     4.0D+00, 2.0D+00, &
     1.0D+00, 5.0D+00, &
    -2.0D+00, 2.0D+00, &
     4.0D+00, 2.0D+00, &
     5.0D+00, 4.0D+00, &
     6.0D+00, 6.0D+00, &
     4.0D+00, 2.0D+00, &
     1.0D+00, 5.0D+00, &
     4.0D+00, 2.0D+00 /), (/ 2, 3, test_num /) )
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE_EDGE_LENGTH_TEST():'
  write ( *, '(a)' ) '  TRIANGLE_EDGE_LENGTH() computes the edge lengths of a triangle.'

  do test = 1, test_num

    t(1:2,1:3) = t_test(1:2,1:3,test)

    call r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

    call triangle_edge_length ( t, edge_length )

    call r8vec_print ( 3, edge_length, '  Edge lengths:' )

  end do

  return
end
subroutine triangle_incircle_test ( )

!*****************************************************************************80
!
!! TRIANGLE_INCIRCLE_TEST() tests TRIANGLE_INCIRCLE();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) pc(2)
  real ( kind = rk8 ) r
  real ( kind = rk8 ), dimension ( 2, 3 ) :: t = reshape ( (/ &
    0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00 /), (/ 2, 3 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE_INCIRCLE_TEST():'
  write ( *, '(a)' ) '  TRIANGLE_INCIRCLE() computes the incircle of a triangle.'

  call r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  call triangle_incircle ( t, r, pc )

  call r8vec_print ( 2, pc, '  Incenter' )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Incircle radius is ', r

  return
end
subroutine triangle_orientation_test ( )

!*****************************************************************************80
!
!! TRIANGLE_ORIENTATION_TEST() tests TRIANGLE_ORIENTATION().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: test_num = 4

  integer triangle_orientation
  integer i
  real ( kind = rk8 ) t(2,3)
  real ( kind = rk8 ), dimension(2,3,test_num) :: t_test = reshape ( (/ &
     4.0D+00,  2.0D+00, &
     1.0D+00,  5.0D+00, &
    -2.0D+00,  2.0D+00, &
     1.0D+00,  5.0D+00, &
     4.0D+00,  2.0D+00, &
     1.0D+00, -1.0D+00, &
     1.0D+00,  5.0D+00, &
     2.0D+00,  7.0D+00, &
     3.0D+00,  9.0D+00, &
     1.0D+00,  5.0D+00, &
     4.0D+00,  2.0D+00, &
     1.0D+00,  5.0D+00 /), (/ 2, 3, test_num /) )
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE_ORIENTATION_TEST():'
  write ( *, '(a)' ) '  TRIANGLE_ORIENTATION() determines the orientation of a triangle.'

  do test = 1, test_num

    t(1:2,1:3) = t_test(1:2,1:3,test)

    i = triangle_orientation ( t )

    call r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

    write ( *, '(a)' ) ''

    if ( i == 0 ) then
      write ( *, '(a)' ) '  The points are counterclockwise.'
    else if ( i == 1 ) then
      write ( *, '(a)' ) '  The points are clockwise.'
    else if ( i == 2 ) then
      write ( *, '(a)' ) '  The points are colinear.'
    else if ( i == 3 ) then
      write ( *, '(a)' ) '  The points are not distinct.'
    else
      write ( *, '(a)' ) '  The return value makes no sense.'
    end if

  end do

  return
end
subroutine triangle_orthocenter_test ( )

!*****************************************************************************80
!
!! TRIANGLE_ORTHOCENTER_TEST() tests TRIANGLE_ORTHOCENTER();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: test_num = 4

  logical flag
  real ( kind = rk8 ) pc(2)
  real ( kind = rk8 ) t(2,3)
  real ( kind = rk8 ), dimension(2,3,test_num) :: t_test = reshape ( (/ &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
     0.0D+00,  1.0D+00, &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
     0.5D+00,  0.86602539D+00, &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
     0.5D+00, 10.0D+00, &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
    10.0D+00,  2.0D+00 /), (/ 2, 3, test_num /) )
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE_ORTHOCENTER_TEST():'
  write ( *, '(a)' ) '  TRIANGLE_ORTHOCENTER() computes the orthocenter of a triangle.'

  do test = 1, test_num

    t(1:2,1:3) = t_test(1:2,1:3,test)

    call r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

    call triangle_orthocenter ( t, pc, flag )

    call r8vec_print ( 2, pc, '  Orthocenter' )

  end do

  return
end
subroutine triangle_point_dist_test ( )

!*****************************************************************************80
!
!! TRIANGLE_POINT_DIST_TEST() tests TRIANGLE_POINT_DIST();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: test_num = 7

  real ( kind = rk8 ) dist
  real ( kind = rk8 ) p(2)
  real ( kind = rk8 ), dimension(2,test_num) :: p_test = reshape ( (/ &
     0.25D+00,   0.25D+00, &
     0.75D+00,   0.25D+00, &
     1.00D+00,   1.00D+00, &
    11.00D+00,   0.50D+00, &
     0.00D+00,   1.00D+00, &
     0.50D+00, -10.00D+00, &
     0.60D+00,   0.60D+00 /), (/ 2, test_num /) )
  real ( kind = rk8 ), dimension ( 2, 3 ) :: t = reshape ( (/ &
    0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00 /), (/ 2, 3 /) )
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE_POINT_DIST_TEST():'
  write ( *, '(a)' ) '  TRIANGLE_POINT_DIST() computes the distance'
  write ( *, '(a)' ) '  from a point to a triangle.'

  call r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           P            DIST'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:2) = p_test(1:2,test)

    call triangle_point_dist ( t, p, dist ) 

    write ( *, '(2x,2f8.3,2x,f8.3)' ) p(1:2), dist

  end do
 
  return
end
subroutine triangle_point_near_test ( )

!*****************************************************************************80
!
!! TRIANGLE_POINT_NEAR_TEST() tests TRIANGLE_POINT_NEAR().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: test_num = 7

  real ( kind = rk8 ) dist
  real ( kind = rk8 ) p(2)
  real ( kind = rk8 ), dimension(2,test_num) :: p_test = reshape ( (/ &
     0.25D+00,   0.25D+00, &
     0.75D+00,   0.25D+00, &
     1.00D+00,   1.00D+00, &
    11.00D+00,   0.50D+00, &
     0.00D+00,   1.00D+00, &
     0.50D+00, -10.00D+00, &
     0.60D+00,   0.60D+00 /), (/ 2, test_num /) )
  real ( kind = rk8 ) pn(2)
  real ( kind = rk8 ), dimension ( 2, 3 ) :: t = reshape ( (/ &
    0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00 /), (/ 2, 3 /) )
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE_POINT_NEAR_TEST():'
  write ( *, '(a)' ) '  TRIANGLE_POINT_NEAR() computes the nearest'
  write ( *, '(a)' ) '  point on a triangle to a given point.'

  call r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           P                PN'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:2) = p_test(1:2,test)

    call triangle_point_near ( t, p, pn, dist ) 

    write ( *, '(2x,2f8.3,2x,2f8.3)' ) p(1:2), pn(1:2)

  end do
 
  return
end
subroutine triangle_quality_test ( )

!*****************************************************************************80
!
!! TRIANGLE_QUALITY_TEST() tests TRIANGLE_QUALITY().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: test_num = 4

  real ( kind = rk8 ) quality
  real ( kind = rk8 ) t(2,3)
  real ( kind = rk8 ), dimension (2,3,test_num) :: t_test = reshape ( (/ &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
     0.0D+00,  1.0D+00, &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
     0.5D+00,  0.86602539D+00, &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
     0.5D+00, 10.0D+00, &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
    10.0D+00,  2.0D+00 /), (/ 2, 3, test_num /) )
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE_QUALITY_TEST():'
  write ( *, '(a)' ) '  TRIANGLE_QUALITY() computes the quality of a triangle.'

  do test = 1, test_num

    t(1:2,1:3) = t_test(1:2,1:3,test)

    call r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

    call triangle_quality ( t, quality )

    write ( *, '(a)' ) ''
    write ( *, '(a,g14.6)' ) '  Quality = ', quality

  end do

  return
end
subroutine triangle_reference_sample_test ( )

!*****************************************************************************80
!
!! TRIANGLE_REFERENCE_SAMPLE_TEST() tests TRIANGLE_REFERENCE_SAMPLE().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: test_num = 10

  real ( kind = rk8 ) p(2)
  real ( kind = rk8 ), dimension(2,3) :: t = reshape ( (/ &
     0.0D+00, 0.0D+00, &
     1.0D+00, 0.0D+00, &
     0.0D+00, 1.0D+00 /), (/ 2, 3 /) )
  integer test
  real ( kind = rk8 ) xsi(3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE_REFERENCE_SAMPLE_TEST():'
  write ( *, '(a)' ) '  TRIANGLE_REFERENCE_SAMPLE() samples the reference triangle.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sample points (X,Y) and (XSI1,XSI2,XSI3) coordinates:'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    call triangle_reference_sample ( 1, p )
    call triangle_xy_to_xsi ( t, p, xsi )
    write ( *, '(2x,2f8.4,4x,3f8.4)' ) p(1:2), xsi(1:3)
  end do

  return
end
subroutine triangle_sample_test ( )

!*****************************************************************************80
!
!! triangle_sample_test() tests triangle_sample().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: test_num = 10

  real ( kind = rk8 ) p(2)
  real ( kind = rk8 ), dimension(2,3) :: t = reshape ( (/ &
     4.0D+00, 2.0D+00, &
     1.0D+00, 5.0D+00, &
    -2.0D+00, 2.0D+00 /), (/ 2, 3 /) )
  integer test
  real ( kind = rk8 ) xsi(3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'triangle_sample_test():'
  write ( *, '(a)' ) '  triangle_sample() samples points from a triangle.'

  call r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sample points (X,Y) and (XSI1,XSI2,XSI3) coordinates:'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    call triangle_sample ( t, 1, p )
    call triangle_xy_to_xsi ( t, p, xsi )
    write ( *, '(2x,2f8.4,4x,3f8.4)' ) p(1:2), xsi(1:3)
  end do

  return
end
subroutine triangle_sample_reflect_test ( )

!*****************************************************************************80
!
!! triangle_sample_reflect_test() tests triangle_sample_reflect().
!
!  Discussion:
!
!    Triangles and point vectors are given as dimension (2,*).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 January 2025
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: n = 1
  integer, parameter :: test_num = 10

  real ( kind = rk8 ) p(2,n)
  real ( kind = rk8 ), dimension(2,3) :: t = reshape ( (/ &
     4.0D+00, 2.0D+00, &
     1.0D+00, 5.0D+00, &
    -2.0D+00, 2.0D+00 /), (/ 2, 3 /) )
  integer test
  real ( kind = rk8 ) xsi(3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'triangle_sample_reflect_test():'
  write ( *, '(a)' ) '  triangle_sample_reflect() samples points from a triangle'
  write ( *, '(a)' ) '  using reflection.'

  call r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sample points (X,Y) and (XSI1,XSI2,XSI3) coordinates:'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    call triangle_sample_reflect ( t, n, p )
    call triangle_xy_to_xsi ( t, p, xsi )
    write ( *, '(2x,2f8.4,4x,3f8.4)' ) p(1:2,1), xsi(1:3)
  end do

  return
end
subroutine triangle_xsi_to_xy_test ( )

!*****************************************************************************80
!
!! TRIANGLE_XSI_TO_XY_TEST() tests TRIANGLE_XSI_TO_XY().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: test_num = 10

  integer i
  real ( kind = rk8 ) p(2)
  real ( kind = rk8 ) p2(2)
  real ( kind = rk8 ), dimension(2,3) :: t = reshape ( (/ &
     4.0D+00, 2.0D+00, &
     1.0D+00, 5.0D+00, &
    -2.0D+00, 2.0D+00 /), (/ 2, 3 /) )
  integer test
  real ( kind = rk8 ) xsi(3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE_XSI_TO_XY_TEST():'
  write ( *, '(a)' ) '  TRIANGLE_XSI_TO_XY() converts XSI to XY coordinates.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We verify that (X,Y) -> (XSI1,XSI2,XSI3) -> (X,Y)'
  write ( *, '(a)' ) '  works properly.'

  call r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sample points:'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    if ( test == 1 ) then
      do i = 1, 2
        p(i) = sum ( t(i,1:3) ) / 3.0D+00
      end do
    else if ( test == 2 ) then
      p(1) = 3.0D+00
      p(2) = 0.0D+00
    else
      call triangle_sample ( t, 1, p )
    end if

    call triangle_xy_to_xsi ( t, p, xsi )

    call triangle_xsi_to_xy ( t, xsi, p2 )

    write ( *, '(2x,2f8.4,4x,3f8.4,2x,2f8.4)' ) p(1:2), xsi(1:3), p2(1:2)

  end do

  return
end
subroutine triangle_xy_to_xsi_test ( )

!*****************************************************************************80
!
!! TRIANGLE_XY_TO_XSI_TEST() tests TRIANGLE_XY_TO_XSI().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: test_num = 10

  integer i
  real ( kind = rk8 ) p(2)
  real ( kind = rk8 ) p2(2)
  real ( kind = rk8 ), dimension(2,3) :: t = reshape ( (/ &
     4.0D+00, 2.0D+00, &
     1.0D+00, 5.0D+00, &
    -2.0D+00, 2.0D+00 /), (/ 2, 3 /) )
  integer test
  real ( kind = rk8 ) xsi(3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE_XY_TO_XSI_TEST():'
  write ( *, '(a)' ) '  TRIANGLE_XY_TO_XSI() converts XY to XSI coordinates.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We verify that (X,Y) -> (XSI1,XSI2,XSI3) -> (X,Y)'
  write ( *, '(a)' ) '  works properly.'

  call r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sample points:'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    if ( test == 1 ) then
      do i = 1, 2
        p(i) = sum ( t(i,1:3) ) / 3.0D+00
      end do
    else if ( test == 2 ) then
      p(1) = 3.0D+00
      p(2) = 0.0D+00
    else
      call triangle_sample ( t, 1, p )
    end if

    call triangle_xy_to_xsi ( t, p, xsi )

    call triangle_xsi_to_xy ( t, xsi, p2 )

    write ( *, '(2x,2f8.4,4x,3f8.4,2x,2f8.4)' ) p(1:2), xsi(1:3), p2(1:2)

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
!  Parameters:
!
!    None
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

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

