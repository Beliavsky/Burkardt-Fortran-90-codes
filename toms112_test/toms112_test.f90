program main

!*****************************************************************************80
!
!! toms112_test() tests toms112().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 November 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'toms112_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test TOMS112().'

  call point_in_polygon_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS112_TEST():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine point_in_polygon_test ( )

!*****************************************************************************80
!
!! POINT_IN_POLYGON_TEST tests POINT_IN_POLYGON.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 November 2016
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
  logical point_in_polygon
  integer test
  real ( kind = rk ) :: x(n) = (/ &
    0.0D+00,  1.0D+00,  2.0D+00,  1.0D+00,  0.0D+00 /)
  real ( kind = rk ) x0
  real ( kind = rk ) :: x0_test(test_num) = (/ &
    1.0D+00, 3.0D+00, 0.0D+00, 0.5D+00 /)
  real ( kind = rk ) :: y(n) = (/ &
    0.0D+00, 0.0D+00, 1.0D+00, 2.0D+00, 2.0D+00 /)
  real ( kind = rk ) y0
  real ( kind = rk ) :: y0_test(test_num) = (/ &
    1.0D+00, 4.0D+00, 2.0D+00, -0.25D+00 /)
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POINT_IN_POLYGON_TEST'
  write ( *, '(a)' ) '  POINT_IN_POLYGON determines if '
  write ( *, '(a)' ) '  a point is in a polygon.'

  call r8vec2_print ( n, x, y, '  The polygon vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          P          Inside'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
 
    x0 = x0_test(test)
    y0 = y0_test(test)

    inside = point_in_polygon ( n, x, y, x0, y0 )

    write ( *, '(2x,2g14.6,4x,l1)' ) x0, y0, inside

  end do
 
  return
end
