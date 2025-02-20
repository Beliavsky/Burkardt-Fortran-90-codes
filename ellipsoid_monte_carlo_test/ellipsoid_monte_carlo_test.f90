program main

!*****************************************************************************80
!
!! ellipsoid_monte_carlo_test() tests ellipsoid_monte_carlo().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 August 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ELLIPSOID_MONTE_CARLO_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test ELLIPSOID_MONTE_CARLO.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ELLIPSOID_MONTE_CARLO_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 uses ELLIPSOID_SAMPLE on a 2D ellipse centered at (0,0).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 August 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 2

  real ( kind = rk ), dimension ( m, m ) :: a = reshape ( (/ &
    9.0, 1.0, &
    1.0, 4.0 /), (/ m, m /) )
  integer e(m)
  integer :: e_test(m,7) = reshape ( (/ &
    0, 0, &
    1, 0, &
    0, 1, &
    2, 0, &
    1, 1, &
    0, 2, &
    3, 0 /), (/ m, 7 /) )
  real ( kind = rk ) ellipsoid_volume
  integer j
  integer n
  real ( kind = rk ), parameter :: r = 2.0D+00
  real ( kind = rk ) result(7)
  real ( kind = rk ), dimension ( m ) :: v = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = rk ), allocatable :: value(:)
  real ( kind = rk ) volume
  real ( kind = rk ), allocatable :: x(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Use ELLIPSOID_SAMPLE to estimate integrals'
  write ( *, '(a)' ) '  in a 2D ellipse x'' * A * x <= r^2.'

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Ellipsoid radius R = ', r
  call r8vec_print ( m, v, '  Ellipsoid center V:' )
  call r8mat_print ( m, m, a, '  Ellipsoid matrix A:' )

  volume = ellipsoid_volume ( m, a, v, r )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Ellipsoid volume = ', volume
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        1              X               Y  ' // &
    '             X^2               XY             Y^2             X^3'
  write ( *, '(a)' ) ' '

  n = 1

  do while ( n <= 65536 )

    allocate ( value(1:n) )
    allocate ( x(1:m,1:n) )

    call ellipsoid_sample ( m, n, a, v, r, x )

    do j = 1, 7

      e(1:m) = e_test(1:m,j)

      call monomial_value ( m, n, e, x, value )

      result(j) = volume * sum ( value(1:n) ) / real ( n, kind = rk )

    end do

    write ( *, '(2x,i8,7(2x,g14.6))' ) n, result(1:7)

    deallocate ( value )
    deallocate ( x )

    n = 2 * n

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 uses ELLIPSOID_SAMPLE on a 2D ellipse centered at (2,3).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 August 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 2

  real ( kind = rk ), dimension ( m, m ) :: a = reshape ( (/ &
    9.0, 1.0, &
    1.0, 4.0 /), (/ m, m /) )
  integer e(m)
  integer :: e_test(m,7) = reshape ( (/ &
    0, 0, &
    1, 0, &
    0, 1, &
    2, 0, &
    1, 1, &
    0, 2, &
    3, 0 /), (/ m, 7 /) )
  real ( kind = rk ) ellipsoid_volume
  integer j
  integer n
  real ( kind = rk ), parameter :: r = 0.5D+00
  real ( kind = rk ) result(7)
  real ( kind = rk ), dimension ( m ) :: v = (/ 2.0D+00, 3.0D+00 /)
  real ( kind = rk ), allocatable :: value(:)
  real ( kind = rk ) volume
  real ( kind = rk ), allocatable :: x(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Use ELLIPSOID_SAMPLE to estimate integrals'
  write ( *, '(a)' ) '  in a 2D ellipse (x-v)'' * A * (x-v) <= r^2.'

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Ellipsoid radius R = ', r
  call r8vec_print ( m, v, '  Ellipsoid center V:' )
  call r8mat_print ( m, m, a, '  Ellipsoid matrix A:' )

  volume = ellipsoid_volume ( m, a, v, r )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Ellipsoid volume = ', volume
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        1              X               Y  ' // &
    '             X^2               XY             Y^2             X^3'
  write ( *, '(a)' ) ' '

  n = 1

  do while ( n <= 65536 )

    allocate ( value(1:n) )
    allocate ( x(1:m,1:n) )

    call ellipsoid_sample ( m, n, a, v, r, x )

    do j = 1, 7

      e(1:m) = e_test(1:m,j)

      call monomial_value ( m, n, e, x, value )

      result(j) = volume * sum ( value(1:n) ) / real ( n, kind = rk )

    end do

    write ( *, '(2x,i8,7(2x,g14.6))' ) n, result(1:7)

    deallocate ( value )
    deallocate ( x )

    n = 2 * n

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 uses ELLIPSOID_SAMPLE on a 3D ellipse centered at (1,2,3).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 August 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 3

  real ( kind = rk ), dimension ( m, m ) :: a = reshape ( (/ &
    9.0, 6.0, 3.0, &
    6.0, 5.0, 4.0, &
    3.0, 4.0, 9.0 /), (/ m, m /) )
  integer e(m)
  integer :: e_test(m,7) = reshape ( (/ &
    0, 0, 0, &
    1, 0, 0, &
    0, 1, 0, &
    0, 0, 1, &
    2, 0, 0, &
    0, 2, 2, &
    0, 0, 3 /), (/ m, 7 /) )
  real ( kind = rk ) ellipsoid_volume
  integer j
  integer n
  real ( kind = rk ), parameter :: r = 0.5D+00
  real ( kind = rk ) result(7)
  real ( kind = rk ), dimension ( m ) :: v = (/ 1.0D+00, 2.0D+00, 3.0D+00 /)
  real ( kind = rk ), allocatable :: value(:)
  real ( kind = rk ) volume
  real ( kind = rk ), allocatable :: x(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Use ELLIPSOID_SAMPLE to estimate integrals'
  write ( *, '(a)' ) '  in a 3D ellipse (x-v)'' * A * (x-v) <= r^2.'

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Ellipsoid radius R = ', r
  call r8vec_print ( m, v, '  Ellipsoid center V:' )
  call r8mat_print ( m, m, a, '  Ellipsoid matrix A:' )

  volume = ellipsoid_volume ( m, a, v, r )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Ellipsoid volume = ', volume
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        1              X               Y  ' // &
    '              Z                X^2            YZ              Z^3'
  write ( *, '(a)' ) ' '

  n = 1

  do while ( n <= 65536 )

    allocate ( value(1:n) )
    allocate ( x(1:m,1:n) )

    call ellipsoid_sample ( m, n, a, v, r, x )

    do j = 1, 7

      e(1:m) = e_test(1:m,j)

      call monomial_value ( m, n, e, x, value )

      result(j) = volume * sum ( value(1:n) ) / real ( n, kind = rk )

    end do

    write ( *, '(2x,i8,7(2x,g14.6))' ) n, result(1:7)

    deallocate ( value )
    deallocate ( x )

    n = 2 * n

  end do

  return
end
