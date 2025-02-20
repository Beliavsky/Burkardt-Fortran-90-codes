program main

!*****************************************************************************80
!
!! spiral_exact_test() tests spiral_exact().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'spiral_exact_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test spiral_exact().'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'spiral_exact_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 generates a field and estimates its range.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  real ( kind = rk ) c
  integer seed
  real ( kind = rk ), allocatable :: u(:)
  real ( kind = rk ), allocatable :: v(:)
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ) xy_hi
  real ( kind = rk ) xy_lo
  real ( kind = rk ), allocatable :: y(:)

  write (  *, '(a)' ) ''
  write (  *, '(a)' ) 'TEST01'
  write (  *, '(a)' ) '  Sample a spiral velocity field and estimate'
  write (  *, '(a)' ) '  the range of the solution values.'

  n = 1000
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  xy_lo = +0.0D+00
  xy_hi = +1.0D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, xy_lo, xy_hi, seed, x )
  call r8vec_uniform_ab ( n, xy_lo, xy_hi, seed, y )
  c = 1.0D+00

  call uv_spiral ( n, x, y, c, u, v )

  write (  *, '(a)' ) ''
  write (  *, '(a)' ) '           Minimum       Maximum'
  write (  *, '(a)' ) ''
  write (  *, '(a,2x,g14.6,2x,g14.6)' ) '  U:  ', minval ( u ), maxval ( u )
  write (  *, '(a,2x,g14.6,2x,g14.6)' ) '  V:  ', minval ( v ), maxval ( v )

  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 generates a field and samples its residuals.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  real ( kind = rk ) c
  real ( kind = rk ), allocatable :: pr(:)
  integer seed
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ) xy_hi
  real ( kind = rk ) xy_lo
  real ( kind = rk ), allocatable :: y(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  Sample a spiral velocity field and estimate the'
  write ( *, '(a)' ) '  range of residuals in the continuity equation.'

  n = 1000
  allocate ( pr(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  xy_lo = +0.0D+00
  xy_hi = +1.0D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, xy_lo, xy_hi, seed, x )
  call r8vec_uniform_ab ( n, xy_lo, xy_hi, seed, y )
  c = 1.0D+00

  call resid_spiral ( n, x, y, c, pr )

  write (  *, '(a)' ) ''
  write (  *, '(a)' ) '           Minimum       Maximum'
  write (  *, '(a)' ) ''
  write (  *, '(a,2x,g14.6,2x,g14.6)' ) '  Pr:  ', minval ( abs ( pr ) ), maxval ( abs ( pr ) )

  deallocate ( pr )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 generates a field on a regular grid and plots it.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: x_num = 21
  integer, parameter :: y_num = 21

  real ( kind = rk ) c
  character ( len = 255 ) header
  integer n
  real ( kind = rk ) s
  real ( kind = rk ), allocatable :: u(:,:)
  real ( kind = rk ), allocatable :: v(:,:)
  real ( kind = rk ), allocatable :: x(:,:)
  real ( kind = rk ) x_hi
  real ( kind = rk ) x_lo
  real ( kind = rk ), allocatable :: y(:,:)
  real ( kind = rk ) y_hi
  real ( kind = rk ) y_lo

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  Generate a spiral velocity field on a regular grid.'
  write ( *, '(a)' ) '  Store in GNUPLOT data and command files.'

  x_lo = 0.0D+00
  x_hi = 1.0D+00

  y_lo = 0.0D+00
  y_hi = 1.0D+00

  allocate ( x(1:x_num,1:y_num) )
  allocate ( y(1:x_num,1:y_num) )

  call grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y )

  n = x_num * y_num
  c = 1.0D+00

  allocate ( u(1:x_num,1:y_num) )
  allocate ( v(1:x_num,1:y_num) )

  call uv_spiral ( n, x, y, c, u, v )

  header = 'spiral'
  s = 0.05D+00
  call spiral_gnuplot ( header, n, x, y, u, v, s )

  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )

  return
end

