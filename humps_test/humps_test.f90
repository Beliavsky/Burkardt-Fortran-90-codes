program main

!*****************************************************************************80
!
!! humps_test() tests humps().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 August 2019
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'humps_test'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test humps.'

  call humps_antideriv_test ( )
  call humps_fun_test ( )
  call humps_deriv_test ( )
  call humps_deriv2_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'humps_test'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine humps_antideriv_test ( )

!*****************************************************************************80
!
!! humps_antideriv_test tests humps_antideriv.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 August 2019
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) humps_antideriv
  integer i
  integer n
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ), allocatable :: y(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'humps_antideriv_test'
  write ( *, '(a)' ) '  Test humps_antideriv.'

  a = 0.0D+00
  b = 2.0D+00
  n = 101
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  do i = 1, n
    x(i) = ( ( n - i ) * a + i * b ) / ( n - 1 )
    y(i) = humps_antideriv ( x(i) )
  end do
  call plot_xy ( n, x, y, 'humps_antideriv' )

  deallocate ( x )
  deallocate ( y )

  return
end
subroutine humps_fun_test ( )

!*****************************************************************************80
!
!! humps_fun_test tests humps_fun.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 August 2019
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) humps_fun
  integer i
  integer n
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ), allocatable :: y(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'humps_fun_test'
  write ( *, '(a)' ) '  Test humps_fun.'

  a = 0.0D+00
  b = 2.0D+00
  n = 101
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  do i = 1, n
    x(i) = ( ( n - i ) * a + i * b ) / ( n - 1 )
    y(i) = humps_fun ( x(i) )
  end do
  call plot_xy ( n, x, y, 'humps_fun' )

  deallocate ( x )
  deallocate ( y )

  return
end
subroutine humps_deriv_test ( )

!*****************************************************************************80
!
!! humps_deriv_test tests humps_deriv.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 August 2019
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) humps_deriv
  integer i
  integer n
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ), allocatable :: y(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'humps_deriv_test'
  write ( *, '(a)' ) '  Test humps_deriv.'

  a = 0.0D+00
  b = 2.0D+00
  n = 101
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  do i = 1, n
    x(i) = ( ( n - i ) * a + i * b ) / ( n - 1 )
    y(i) = humps_deriv ( x(i) )
  end do
  call plot_xy ( n, x, y, 'humps_deriv' )

  deallocate ( x )
  deallocate ( y )

  return
end
subroutine humps_deriv2_test ( )

!*****************************************************************************80
!
!! humps_deriv2_test tests humps_deriv2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 August 2019
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) humps_deriv2
  integer i
  integer n
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ), allocatable :: y(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'humps_deriv2_test'
  write ( *, '(a)' ) '  Test humps_deriv2.'

  a = 0.0D+00
  b = 2.0D+00
  n = 101
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  do i = 1, n
    x(i) = ( ( n - i ) * a + i * b ) / ( n - 1 )
    y(i) = humps_deriv2 ( x(i) )
  end do
  call plot_xy ( n, x, y, 'humps_deriv2' )

  deallocate ( x )
  deallocate ( y )

  return
end

