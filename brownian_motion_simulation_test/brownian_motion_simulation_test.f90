program main

!*****************************************************************************80
!
!! brownian_motion_simulation_test() tests brownian_motion_simulation().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) d
  real ( kind = rk ), allocatable :: dsq(:,:)
  character ( len = 255 ) header
  integer k
  integer m
  integer n
  real ( kind = rk ) t
  real ( kind = rk ), allocatable :: x(:,:)

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'brownian_motion_simulation_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test brownian_motion_simulation().'
!
!  Compute the path of a particle undergoing Brownian motion.
!
  do m = 1, 2

    n = 1001
    d = 10.0D+00
    t = 1.0D+00
    allocate ( x(1:m,1:n) )
    call brownian_motion_simulation ( m, n, d, t, x )
    if ( m == 1 ) then
      header = 'motion_1d'
    else if ( m == 2 ) then
      header = 'motion_2d'
    end if
    call brownian_motion_display ( m, n, x, header )
    deallocate ( x )

  end do
!
!  Estimate the average displacement of the particle from the origin
!  as a function of time.
!
  do m = 1, 3

    k = 40
    n = 1001
    d = 10.0D+00
    t = 1.0D+00

    allocate ( dsq(1:k,1:n) )
    call brownian_displacement_simulation ( k, n, m, d, t, dsq )
    if ( m == 1 ) then
      header = 'displacement_1d'
    else if ( m == 2 ) then
      header = 'displacement_2d'
    else if ( m == 3 ) then
      header = 'displacement_3d'
    end if
    call brownian_displacement_display ( k, n, m, d, t, dsq, header )
    deallocate ( dsq )

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'brownian_motion_simulation_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( );

  stop 0
end
