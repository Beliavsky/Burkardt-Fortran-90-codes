program main

!*****************************************************************************80
!
!! ball_grid_test() tests ball_grid().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ball_grid_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the BALL_GRID library.'

  call ball_grid_test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BALL_GRID_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine ball_grid_test01 ( )

!*****************************************************************************80
!
!! ball_grid_test01() tests ball_grid().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: bg(:,:)
  real ( kind = rk ) c(3)
  character ( len = 80 ) filename
  integer n
  integer ng
  real ( kind = rk ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  BALL_GRID can define a grid of points'
  write ( *, '(a)' ) '  with N+1 points on a horizontal or vertical radius,'
  write ( *, '(a)' ) '  based on any ball.'

  n = 10
  r = 2.0D+00
  c(1) = 1.0D+00
  c(2) = 5.0D+00
  c(3) = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  We use N = ', n
  write ( *, '(a,g14.6)' ) '  Radius R = ', r
  write ( *, '(a,g14.6,a,g14.6,a,g14.6,a)' ) &
    '  Center C = (', c(1), ',', c(2), ',', c(3), ')'

  call ball_grid_count ( n, r, c, ng )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of grid points will be ', ng

  allocate ( bg(1:3,1:ng) )

  call ball_grid ( n, r, c, ng, bg )

  call r83vec_print_part ( ng, bg, 20, '  Part of the grid point array:' )

  filename = 'ball_grid_test01.xyz'

  call r8mat_write ( filename, 3, ng, bg )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to the file "' // trim ( filename ) // '".'

  deallocate ( bg )

  return
end
 
