program main

!*****************************************************************************80
!
!! PYRAMID_GRID_TEST tests the PYRAMID_GRID library.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 August 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PYRAMID_GRID_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the PYRAMID_GRID library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PYRAMID_GRID_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests PYRAMID_GRID_SIZE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 August 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer ng
  integer pyramid_grid_size

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  PYRAMID_GRID_SIZE determines the size of a'
  write ( *, '(a)' ) '  pyramid grid with N+1 points along each edge.'

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   N    Size'
  write ( *, '(a)' ) ''
  do n = 0, 10
    ng = pyramid_grid_size ( n )
    write ( *, '(2x,i2,2x,i6)' ) n, ng
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests PYRAMID_UNIT_GRID.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 August 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer ng
  real ( kind = rk ), allocatable :: pg(:,:)
  integer pyramid_grid_size

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  PYRAMID_UNIT_GRID determines a unit pyramid'
  write ( *, '(a)' ) '  grid with N+1 points along each edge.'

  n = 4
  call r8_print ( n, '  Grid parameter N:' )

  ng = pyramid_grid_size ( n )
  call r8_print ( ng, '  Grid size NG:' )

  allocate ( pg(1:3,1:ng) )

  call pyramid_unit_grid ( n, ng, pg )

  call r8mat_transpose_print ( 3, ng, pg, '  Pyramid grid points:' )

  deallocate ( pg )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests PYRAMID_UNIT_GRID_PLOT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 August 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 255 ) header
  integer n
  integer ng
  real ( kind = rk ), allocatable :: pg(:,:)
  integer pyramid_grid_size

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  PYRAMID_UNIT_GRID_PLOT plots a unit pyramid'
  write ( *, '(a)' ) '  grid with N+1 points along each edge.'

  n = 5
  call r8_print ( n, '  Grid parameter N:' )

  ng = pyramid_grid_size ( n )
  call r8_print ( ng, '  Grid size NG:' )

  allocate ( pg(1:3,1:ng) )

  call pyramid_unit_grid ( n, ng, pg )

  header = 'pyramid_unit'
  call pyramid_unit_grid_plot ( n, ng, pg, header )

  deallocate ( pg )

  return
end
