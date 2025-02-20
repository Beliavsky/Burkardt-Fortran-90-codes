program main

!*****************************************************************************80
!
!! SPHERE_FIBONACCI_GRID_TEST tests the SPHERE_FIBONACCI_GRID library.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 May 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'SPHERE_FIBONACCI_GRID_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SPHERE_FIBONACCI_GRID library.'

  call sphere_fibonacci_grid_points_test ( )
  call sphere_fibonacci_grid_display_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'SPHERE_FIBONACCI_GRID_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine sphere_fibonacci_grid_points_test ( )

!*****************************************************************************80
!
!! SPHERE_FIBONACCI_GRID_POINTS_TEST tests SPHERE_FIBONACCI_GRID_POINTS.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 May 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 255 ) filename
  integer ng
  real ( kind = rk ), allocatable :: xg(:,:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'SPHERE_FIBONACCI_GRID_POINTS_TEST'
  write ( *, '(a)' ) '  SPHERE_FIBONACCI_GRID_POINTS computes points on a sphere'
  write ( *, '(a)' ) '  that lie on a Fibonacci spiral.'

  ng = 1000
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Number of points NG = ', ng

  allocate ( xg(1:3,1:ng) )

  call sphere_fibonacci_grid_points ( ng, xg )

  call r8mat_transpose_print_some ( 3, ng, xg, 1, 1, 3, 10, &
    '  Part of the grid array:' )
!
!  Write the nodes to a file.
!
  filename = 'sphere_fibonacci_grid_n1000.xyz'

  call r8mat_write ( filename, 3, ng, xg )

  deallocate ( xg )

  return
end
subroutine sphere_fibonacci_grid_display_test ( )

!*****************************************************************************80
!
!! SPHERE_FIBONACCI_GRID_DISPLAY_TEST tests SPHERE_FIBONACCI_GRID_DISPLAY.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 May 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ng
  character ( len = 255 ) prefix
  real ( kind = rk ), allocatable :: xg(:,:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'SPHERE_FIBONACCI_GRID_DISPLAY_TEST'
  write ( *, '(a)' ) '  SPHERE_FIBONACCI_GRID_DISPLAY displays points'
  write ( *, '(a)' ) '  on a sphere that lie on a Fibonacci spiral.'

  ng = 1000
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Number of points NG = ', ng

  allocate ( xg(1:3,1:ng) )

  call sphere_fibonacci_grid_points ( ng, xg )
!
!  Display the nodes on a sphere.
!
  prefix = 'sphere_fibonacci_grid_n1000'

  call sphere_fibonacci_grid_display ( ng, xg, prefix )

  deallocate ( xg )

  return
end
